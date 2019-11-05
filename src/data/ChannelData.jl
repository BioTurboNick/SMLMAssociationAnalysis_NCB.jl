"""
    ChannelData

Stores the results from a single channel.
"""
struct ChannelData
    channelname::AbstractString
    molecules::Vector{<:DataEntity}
    neighbormolecules::Vector{<:DataEntity}
end

"""
    ChannelDataOptimizing

Stores the results from many molecule lists. Intended for parameter optimizing only.
"""
struct ChannelDataOptimizing
    channelname::AbstractString
    molecules::Vector{Vector{<:DataEntity}}
    neighbormolecules::Vector{<:DataEntity}
end
