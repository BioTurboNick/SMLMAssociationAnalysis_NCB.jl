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
export ChannelDataOptimizing, ResultOptimizing
export getmolecules, exclusivenearestneighbors, montecarloaffinity
export moleculesplot, localizationsplot, insetplot, distanceprobabilityplot
export localizationsplot_forprint
export moleculesinsetplot, localizationsinsetplot
export moleculesinsetplot_forprint, localizationsinsetplot_forprint
end
