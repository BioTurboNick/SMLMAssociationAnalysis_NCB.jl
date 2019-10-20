using Profile

"""
    groupby_dbscan_temporallimit()

Group nearby localizations separated by time into `Molecule`s

Identifies local maxima of `Localization`s temporally distant from other points by `t_off` frames. The effect is to
consider flashes correlated in time as originating from the same molecule. Only localizations within `radius` (same
units as position) of the cluster center are joined into the molecule. The algorithm is repeated on the remaining
unjoined localizations.
"""
function groupby_localmax_temporallimit(localizations::Vector{Localization},
    radius, t_off)

    molecules = Molecule[]


    println("groupby")

    while length(localizations) > 0
        println("iteration")
        @time neighborsdict = findtemporalneighbors(localizations, radius, t_off) # very significant allocations but fast

        @time localmaxima, localmaxima_map = findlocalmaxima(localizations, neighborsdict) # significant allocations and only a little slow

        @time setdiff!(localizations, localmaxima) # remaining = filter(l -> l ∉ localmaxima, localizations) # slooooooooow step 1240 seconds!!!! But virtually no allocations
        #Profile.print(noisefloor = 2.0, mincount = 100)
        #Profile.print(format = :flat, sortedby = :count, noisefloor = 2.0, mincount = 1000)
        @time remaining_localmaxima_map = IdDict([Pair(k, localmaxima_map[k]) for k ∈ setdiff(keys(localmaxima_map), localmaxima)]) # filter(kv -> kv.first ∉ localmaxima, localmaxima_map) # sloooooow 820 seconds!!! 350k allocations, but only 19 MB

        @time newmolecules = map(l -> Molecule(l), localmaxima)  # fast, 140k allocations and 9 MB

        @time buildmolecules!(newmolecules, localizations, remaining_localmaxima_map, radius) #buildmolecules!(newmolecules, remaining, remaining_localmaxima_map, radius) # slooow 472 seconds, 100M allocations! 100 GB!

        #@time localizations = remaining    # instant
        @time append!(molecules, newmolecules) #instant
    end

    molecules
end

"""
    merge_close_molecules

Merge closely clustered molecules (which may be physically associated) into single molecules.
"""
function merge_close_molecules(molecules::Vector{Molecule}, radius)
    mergedmolecules = Molecule[]

    while length(molecules) > 0
        neighborsdict = findneighbors(molecules, radius)

        localmaxima, localmaxima_map = findlocalmaxima(molecules, neighborsdict)

        #remaining = filter(l -> l ∉ localmaxima, molecules)
        setdiff!(molecules, localmaxima)
        #remaining_localmaxima_map = filter(kv -> kv.first ∉ localmaxima, localmaxima_map)
        remaining_localmaxima_map = IdDict([Pair(k, localmaxima_map[k]) for k ∈ setdiff(keys(localmaxima_map), localmaxima)])

        buildmolecules!(localmaxima, molecules, remaining_localmaxima_map, radius)

        #molecules = remaining
        append!(mergedmolecules, localmaxima)
    end

    mergedmolecules
end

"""
    findtemporalneighbors()

Find `Localiation`s that occur near each other within `t_off` frames of each other.
"""
function findtemporalneighbors(localizations::Vector{Localization}, radius, t_off)
        neighborsdict = findneighbors(localizations, radius)

        function istemporallyclose(localization1, localization2)
            if localization1.frame < localization2.frame
                loc1end = localization1.frame + localization2.length - 1
                return localization2.frame - loc1end ≤ t_off
            else
                loc2end = localization2.frame + localization2.length - 1
                return localization1.frame - loc2end ≤ t_off
            end
        end

        foreach(l -> filter!(n -> istemporallyclose(n, l), neighborsdict[l]), localizations)
        neighborsdict
end

"""
    findneighbors()

Find `Localiation`s or `Molecule`s within `radius` of each other.
"""
function findneighbors(points::Vector{<:DataEntity}, radius)
    coordinates = extractcoordinates(points)

    neighbortree = BallTree(coordinates)
    ineighbors = inrange(neighbortree, coordinates, radius, true)
    neighbors = map(i -> points[i], ineighbors)
    IdDict(zip(points, neighbors))
end

findnearestneighbor(neighbors::Vector{<:DataEntity}, molecule::Molecule) =
    findnearestneighbor(neighbors, molecule.group)

"""
    findnearestneighbor()

Find the closest `Localiation` or `Molecule` to another.
"""
function findnearestneighbor(neighbors::Vector{<:DataEntity}, point::DataEntity)
    coordinates = extractcoordinates(point)
    neighborcoordinates = extractcoordinates(neighbors)

    neighbortree = BallTree(neighborcoordinates)
    inearestneighbor, nearestneighbor_distance = nn(neighbortree, coordinates)

    neighbors[inearestneighbor], nearestneighbor_distance
end

"""
    findnearestunvisitedmaxneighbor()

Find the nearest neighbor that hadn't been previously visited.
"""
function findnearestunvisitedmaxneighbor(entity, neighborsdict, visited)
    neighbor_counts = length.(values(neighborsdict))
    maxneighbor_indexes = argmaxall(neighbor_counts)
    maxneighbors = collect(keys(neighborsdict))[maxneighbor_indexes]

    # if there are multiple neighbors with max count, choose closest
    if length(maxneighbors) > 1
        unvisited_maxneighbors = setdiff(maxneighbors, visited)

        length(unvisited_maxneighbors) > 0 &&
            return first(findnearestneighbor(unvisited_maxneighbors, entity))

        return entity
    end

    return first(maxneighbors)
end

"""
    create_neighborsofneighborsdict()

Create a dictionary containing the 2nd-level neighbors of each `Localization` or `Molecule`.
"""
function create_neighborsofneighborsdict(entities::Vector{T}, neighborsdict::IdDict{T, Vector{T}}) where T <: DataEntity
    neighborsofneighborsdict = IdDict{T, IdDict{T, Vector{T}}}()

    foreach(l -> neighborsofneighborsdict[l] = IdDict{T, Vector{T}}(), entities)

    for neighborlistpair ∈ neighborsdict
        entity = neighborlistpair.first
        neighbors = neighborlistpair.second
        foreach(n -> neighborsofneighborsdict[n][entity] = neighbors, neighbors)
    end

    neighborsofneighborsdict
end

"""
    findlocalmaxima()

Find the `Localization`s or `Molecule`s corresponding to local density maxima.
"""
function findlocalmaxima(entities::Vector{T}, neighborsdict::IdDict{T, Vector{T}}) where T <: DataEntity
    localmaxima = Vector{T}()

    # map of each localization to the index of its associated local maximum
    localmaxima_map = IdDict{T, Int}()

    updatemap(chain, i) =
        foreach(x -> localmaxima_map[x] = i, chain)

    isvisited(entity) =
        haskey(localmaxima_map, entity)

    neighborsofneighborsdict = create_neighborsofneighborsdict(entities, neighborsdict)

    for entity ∈ entities
        isvisited(entity) && continue

        visited_chain = Vector{T}()

        while true
            push!(visited_chain, entity)

            localneighborsdict = neighborsofneighborsdict[entity]

            maxneighbor = findnearestunvisitedmaxneighbor(entity, localneighborsdict, visited_chain)

            if (reached_maximum = maxneighbor === entity)
                push!(localmaxima, entity)
                updatemap(visited_chain, entity.index)
                break
            else
                if (joined_chain = isvisited(maxneighbor))
                    localmaxindex = localmaxima_map[maxneighbor]
                    push!(localmaxima_map, maxneighbor => localmaxindex)
                    updatemap(visited_chain, localmaxindex)
                    break
                end

                entity = maxneighbor
            end
        end
    end
    (localmaxima, localmaxima_map)
end

# Is there anywhere that a comprehension would be good?

"""
    buildmolecules!

Build grouped molecules from the molecule seeds and the other entities (`Localization`s or `Molecule`s) to join.
"""
function buildmolecules!(seedmolecules::Vector{Molecule}, entities::Vector{T},
    localmaxima_map::IdDict{T, Int}, radius) where T <: DataEntity

    moleculeneighbors_map = invertdictionary(localmaxima_map)

    for molecule ∈ seedmolecules
        haskey(moleculeneighbors_map, molecule.index) || continue

        clusterneighbors = moleculeneighbors_map[molecule.index]

        while length(clusterneighbors) > 0
            nearestneighbor, distance = findnearestneighbor(clusterneighbors, molecule)

            isinrange = distance < radius + nearestneighbor.accuracy + molecule.accuracy

            isinrange || break

            push!(molecule, nearestneighbor)
            deleteat!(clusterneighbors, findfirst(x -> x === nearestneighbor, clusterneighbors))
            deleteat!(entities, findfirst(x -> x === nearestneighbor, entities))
        end
    end

    seedmolecules
end
