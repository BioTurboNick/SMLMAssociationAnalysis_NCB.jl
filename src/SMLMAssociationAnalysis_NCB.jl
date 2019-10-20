module SMLMAssociationAnalysis_NCB

using Distributed
using LocalizationMicroscopy
using NearestNeighbors
using Plots
using SharedArrays
using Statistics
using Tau

include("data/Molecule.jl")
include("data/ChannelData.jl")
include("data/Result.jl")
include("utilities.jl")
include("analysis/process.jl")
include("analysis/grouping.jl")
include("analysis/association.jl")
include("plots.jl")

export ChannelData, Molecule, Result
export getmolecules, exclusivenearestneighbors, montecarloaffinity

end
