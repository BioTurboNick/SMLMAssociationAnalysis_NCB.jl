"""
    exclusivenearestneighbors

Exclusively pair points between two provided lists.

Finds the closest of all pairs, removes both, and iterates on the remainder until all pairs have been created.

NOTE: Assumes that your data is far away from the type's maximum value.
"""
function exclusivenearestneighbors1(points1::Vector{T}, points2::Vector{T}) where T <: DataEntity
    points2larger = length(points2) > length(points1)

    (points1, points2) = points2larger ?
        (points2, copy(points1)) :
        (points1, copy(points2))

    coordinates1 = extractcoordinates(points1)
    coordinates1 = ndims(coordinates1) == 1 ?
        make2d(coordinates1) :
        coordinates1

    neighbors1 = T[]
    neighbors2 = T[]
    distances = eltype(coordinates1)[]

    neighbortree = BallTree(coordinates1)
    coordinates1used = falses(length(points1))

    while !isempty(points2)
        coordinates2 = extractcoordinates(points2)
        coordinates2used = falses(length(points2))
        inearestneighbors, nearestneighbor_distances = nn(neighbortree, coordinates2, i -> coordinates1used[i])
        while any(nearestneighbor_distances .< typemax(eltype(coordinates1)))
            mindist, mini = findmin(nearestneighbor_distances)
            nearestneighbor_distances[mini] = typemax(eltype(coordinates1))
            coordinates2used[mini] && continue
            inearestneighbor = inearestneighbors[mini]
            coordinates1used[inearestneighbor] && continue
            coordinates2used[mini] = coordinates1used[inearestneighbor] = true
            push!(neighbors1, points1[inearestneighbor])
            push!(neighbors2, points2[mini])
            push!(distances, mindist)
        end
        deleteat!(points2, coordinates2used)
    end

    if points2larger
        neighbors2, neighbors1, distances
    else
        neighbors1, neighbors2, distances
    end
end

function exclusivenearestneighbors(points1::Vector{T}, points2::Vector{T}) where T <: DataEntity
    points2larger = length(points2) > length(points1)

    if points2larger
        points1, points2 = points2, points1
    end

    coordinates1 = extractcoordinates(points1)

    if ndims(coordinates1) == 1
        coordinates1 = make2d(coordinates1) # ensures array is 2d, currently needed for Tree
    end

    neighbors = Vector{T}()
    distances = Vector{eltype(coordinates1)}()

    for point ∈ points2
        coordinates = extractcoordinates(point)

        neighbortree = KDTree(coordinates1)
        inearestneighbor, nearestneighbor_distance = nn(neighbortree, coordinates)

        neighbor = points1[inearestneighbor]

        coordinates1[:,inearestneighbor] .= typemax(eltype(coordinates1))
        push!(neighbors, neighbor)
        push!(distances, nearestneighbor_distance)
    end

    if points2larger
        points2, neighbors, distances
    else
        neighbors, points2, distances
    end
end

"""
    montecarloaffinity()

Evaluate the probability of chance association across all molecules by randomizing their neighbors.
"""
function montecarloaffinity(molecules1::Vector{T}, molecules2::Vector{T}, ch1_neighbors::Vector{T},
                            ch2_neighbors::Vector{T}, distances, maxbindingdistance, rangefactor) where T <: DataEntity
    coordinates1, coordinates2 = extractcoordinates.([molecules1, molecules2])

    percentileranks = ones(length(ch1_neighbors))

    neighborcoordinates1, neighborcoordinates2 = extractcoordinates.([ch1_neighbors, ch2_neighbors])

    centerxcoordinates = mean([neighborcoordinates1[1, :]'; neighborcoordinates2[1, :]'], dims = 1)
    centerycoordinates = mean([neighborcoordinates1[2, :]'; neighborcoordinates2[2, :]'], dims = 1)
    centerzcoordinates = mean([neighborcoordinates1[3, :]'; neighborcoordinates2[3, :]'], dims = 1)
    centercoordinates = [centerxcoordinates; centerycoordinates; centerzcoordinates]

    localradius = rangefactor * maxbindingdistance

    neighbor1tree, neighbor2tree = BallTree.([coordinates1, coordinates2])
    ineighbors1, ineighbors2 = inrange.([neighbor1tree, neighbor2tree], Ref(centercoordinates), localradius, true)
    nlocalmolecules1, nlocalmolecules2 = [length.(x) for x ∈ [ineighbors1, ineighbors2]]

    percentileranks = pmap(localmontecarlo, nlocalmolecules1, nlocalmolecules2, distances, (localradius for i = 1:length(ch1_neighbors)), (10000 for i = 1:length(ch1_neighbors)), on_error = identity)

    return percentileranks
end # since number might be finite, could potentially generate say 100,000 permutations once and use for all with same counts...

"""
    localmontecarlo()

Evaluate the probabilty of chance association of a molecule pair by randomizing its neighbors.
"""
function localmontecarlo(nlocalmolecules1, nlocalmolecules2, testdistance, radius, iterations) where T <: DataEntity
    (nlocalmolecules1 == 0 || nlocalmolecules2 == 0) && return 1.0
    (testdistance > 2 * radius) && return 1.0

    randomcoordinates1 = [randomcoordinates2d(nlocalmolecules1, radius) for i ∈ 1:iterations]
    randomcoordinates2 = [randomcoordinates2d(nlocalmolecules2, radius) for i ∈ 1:iterations]
    randomtrees = nlocalmolecules1 > nlocalmolecules2 ? KDTree.(randomcoordinates1) : KDTree.(randomcoordinates2)
    randomcoordinates = nlocalmolecules1 > nlocalmolecules2 ? randomcoordinates2 : randomcoordinates1
    mindistance = (nn(randomtrees[i], randomcoordinates[i]) |> last for i ∈ 1:iterations) .|> minimum
    percentilerank = count(mindistance .≤ testdistance) / iterations

    return percentilerank
end
