"""
    ChannelData

Stores the results from a single channel.
"""
struct ChannelData
    channelname::AbstractString
    molecules::Vector{<:DataEntity}
    neighbormolecules::Vector{<:DataEntity}
end
