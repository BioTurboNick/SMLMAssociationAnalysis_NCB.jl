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

"""
    ResultOptimizing

Stores the results from an experiment with extra data for optimization purposes.
"""
struct ResultOptimizing
    projectdirname::AbstractString
    experimentdirname::AbstractString
    replicate::Int
    samplename::AbstractString
    cell::Int
    channels::Vector{ChannelDataOptimizing}
    distances::Vector{AbstractFloat}
    mediandistance::AbstractFloat
    percentileranks::Vector{Float64}
    percentileranks_by_localdensity::Vector{Vector{Float64}}
end