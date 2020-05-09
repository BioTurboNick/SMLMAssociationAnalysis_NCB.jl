"""
    Result

Stores the results from an experiment.
"""
struct NewResult
    projectdirname::AbstractString
    experimentdirname::AbstractString
    replicate::Int
    samplename::AbstractString
    cell::Int
    channels::Vector{ChannelData}
    distances::Vector{AbstractFloat}
    mediandistance::AbstractFloat
    percentileranks::Vector{AbstractFloat}
    positivecontrol_percentileranks::Vector{AbstractFloat} # stores the percentile ranks obtained when simulating 100% bound from this cell
    negativecontrol_percentileranks::Vector{AbstractFloat} # stores the percentile ranks obtained when simulating 0% bound from this cell
end

struct Result
    projectdirname::AbstractString
    experimentdirname::AbstractString
    replicate::Int
    samplename::AbstractString
    cell::Int
    channels::Vector{ChannelData}
    distances::Vector{AbstractFloat}
    mediandistance::AbstractFloat
    percentileranks::Vector{AbstractFloat}
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
    percentileranks_by_localdensity::Vector{Vector{Float64}}
end
