"""
    argmaxall()

Return the index of the maximum element in a collection. If there are multiple, all tied indexes are returned.
"""
function argmaxall(vals)
    maxvalue = maximum(vals)
    indexedvals = collect(zip(vals, eachindex(vals)))
    maxindexedvals = filter(v -> v[1] == maxvalue, indexedvals)
    last.(maxindexedvals)
end

"""
    invertdictionary()

Create a new dictionary using the values of an existing dictionary as keys and the values of the existing dictionary as
values.
"""
function invertdictionary(dictionary::IdDict{S, T}) where {S, T}
    newkeys = unique(values(dictionary))
    inverted = Dict{T, Vector{S}}()
    foreach(k -> inverted[k] = Vector{S}(), newkeys)
    foreach(kv -> push!(inverted[kv.second], kv.first), dictionary)
    inverted
end

"""
    inrange()

Evaluate if a value lies within the range provided, inclusive.
"""
inrange(value, lowerlimit, upperlimit) =
    value ≥ lowerlimit && value ≤ upperlimit

"""
    randomcoordinates()

Gernate random coordinates within a radius.
"""
function randomcoordinates2d(count, radius)
    angles = rand(count) * τ
    radii = sqrt.(rand(count) .* (radius ^ 2)) # required for a uniform distribution around the circle
    x = radii .* cos.(angles)
    y = radii .* sin.(angles)
    z = zeros(count)
    return [x'; y'; z']
end
