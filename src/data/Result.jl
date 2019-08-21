"""
    Result

Stores the results from an experiment.
"""
struct Result
    projectdirname::AbstractString
    experimentdirname::AbstractString
    replicate::Int
    samplename::AbstractString
    cell::Int
    channels::Vector{ChannelData}
    distances::Vector{AbstractFloat}
    mediandistance::AbstractFloat
    percentileranks::Vector{Float64}
end
