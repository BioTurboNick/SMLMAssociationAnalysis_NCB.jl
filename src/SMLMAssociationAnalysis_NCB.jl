module SMLMAssociationAnalysis_NCB

using JLD2
using FileIO
using LocalizationMicroscopy
using NearestNeighbors
using Plots
using Printf
using SimpleANOVA
using Statistics
using StatsPlots

include(raw"data\Molecule.jl")
include(raw"data\ChannelData.jl")
include(raw"data\Result.jl")
include("utilities.jl")
include(raw"analysis\process.jl")
include(raw"analysis\grouping.jl")
include(raw"analysis\association.jl")
include("plots.jl")

export ChannelData, Molecule, Result
export getmolecules, exclusivenearestneighbors, montecarloaffinity

end
