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
    locs = copy(localizations)

    while length(locs) > 0
        neighborsdict = findtemporalneighbors(locs, radius, t_off) # 9 seconds, 5 million allocations, 5.6 GB, 35.27% GC time
        localmaxima, localmaxima_map = findlocalmaxima(locs, neighborsdict) # 19 seconds, 13 million allocations, 2 GB, 15.90% GC time
        setdiff!(locs, localmaxima)
        remaining_localmaxima_map = IdDict(Pair(k, localmaxima_map[k]) for k ∈ setdiff(keys(localmaxima_map), localmaxima)) # 1 second 350 thousand allocations, 32 MB
        newmolecules = map(l -> Molecule(l), localmaxima)
        buildmolecules!(newmolecules, locs, remaining_localmaxima_map, radius) # 120 seconds, 3.5 million allocations, 1 GB, 0.2% GC time
        append!(molecules, newmolecules)
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
        setdiff!(molecules, localmaxima)
        remaining_localmaxima_map = IdDict([Pair(k, localmaxima_map[k]) for k ∈ setdiff(keys(localmaxima_map), localmaxima)])
        buildmolecules!(localmaxima, molecules, remaining_localmaxima_map, radius)
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
        setdiff!(maxneighbors, visited)
        return length(maxneighbors) > 0 ?
            first(findnearestneighbor(maxneighbors, entity)) :
            entity
    end

    return first(maxneighbors)
end

"""
    create_neighborsofneighborsdict()

Create a dictionary containing the 2nd-level neighbors of each `Localization` or `Molecule`.
"""
function create_neighborsofneighborsdict(entities::Vector{T}, neighborsdict::IdDict{T, Vector{T}}) where T <: DataEntity
    neighborsofneighborsdict = IdDict{T, IdDict{T, Vector{T}}}()
    foreach(e -> neighborsofneighborsdict[e] = IdDict(Pair(k, neighborsdict[k]) for k ∈ neighborsdict[e]), entities)
    neighborsofneighborsdict
end

"""
    findlocalmaxima()

Find the `Localization`s or `Molecule`s corresponding to local density maxima.
"""
function findlocalmaxima(entities::Vector{T}, neighborsdict::IdDict{T, Vector{T}}) where T <: DataEntity
    localmaxima = T[]

    # map of each localization to the index of its associated local maximum
    localmaxima_map = IdDict{T, Int}()

    updatemap(chain, i) =
        foreach(x -> localmaxima_map[x] = i, chain)

    isvisited(entity) =
        haskey(localmaxima_map, entity)

    neighborsofneighborsdict = create_neighborsofneighborsdict(entities, neighborsdict)

    for entity ∈ entities
        isvisited(entity) && continue

        visited_chain = T[]

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

"""
    buildmolecules!

Build grouped molecules from the molecule seeds and the other entities (`Localization`s or `Molecule`s) to join.
"""
function buildmolecules!(seedmolecules::Vector{Molecule}, entities::Vector{T},
    localmaxima_map::IdDict{T, Int}, radius) where T <: DataEntity

    moleculeneighbors_map = invertdictionary(localmaxima_map)
    entityindexes = IdDict(zip(entities, eachindex(entities)))

    todelete = Int[]

    for molecule ∈ seedmolecules
        haskey(moleculeneighbors_map, molecule.index) || continue

        clusterneighbors = moleculeneighbors_map[molecule.index]
        clusterneighborcoordinates = extractcoordinates(clusterneighbors)
        clusterneighborsused = falses(length(clusterneighbors))
        clusterneighbortree = BallTree(clusterneighborcoordinates)
        coordinates = extractcoordinates(molecule)

        for i ∈ eachindex(clusterneighbors)
            inearestneighbor, nearestneighbor_distance = nn(clusterneighbortree, coordinates, j -> clusterneighborsused[j])
            nearestneighbor = clusterneighbors[inearestneighbor]

            isinrange = nearestneighbor_distance < radius + nearestneighbor.accuracy + molecule.accuracy

            isinrange || break

            push!(molecule, nearestneighbor)
            clusterneighborsused[inearestneighbor] = true
            push!(todelete, entityindexes[nearestneighbor])
            # try batching this one more time; filter not efficient
        end

    end

    sort!(todelete)
    deleteat!(entities, todelete)

    seedmolecules
end
