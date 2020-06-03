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
    percentileranks::Vector{AbstractFloat}
    positivecontrol_distances::Vector{AbstractFloat} # stores the distances obtained when simulating 100% bound from this cell
    negativecontrol_distances::Vector{Vector{AbstractFloat}} # stores the distances obtained when simulating 0% bound from this cell
    positivecontrol_percentileranks::Vector{AbstractFloat} # stores the percentile ranks obtained when simulating 100% bound from this cell
    negativecontrol_percentileranks::Vector{Vector{AbstractFloat}} # stores the percentile ranks obtained when simulating 0% bound from this cell
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

"""
    ResultSimulate

Stores the results from an experiment, for the simulation method (without positive/negative controls)
"""
struct ResultSimulate
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
