module SMLMAssociationAnalysis_NCB

using LocalizationMicroscopy
using NearestNeighbors
using Plots

include(raw"data\Molecule.jl")
include(raw"data\ChannelData.jl")
include(raw"data\Result.jl")
include("utilities.jl")
include(raw"analysis\process.jl")
include(raw"analysis\grouping.jl")
include(raw"analysis\association.jl")
include("plots.jl")

end
