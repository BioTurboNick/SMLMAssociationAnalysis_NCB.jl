#=
Remove localizations from the molecules to reduce results file size.
=#

using FileIO
using SMLMAssociationAnalysis_NCB.jl

outputdir = "../../output"
datapath = joinpath(outputdir, "resultspos.jld2")

experimentresults = load(datapath)["experimentresults"]

experimentresults1 = [[[Result(r, [ChannelData(c.channelname, ResultMolecule.(c.molecules), ResultMolecule.(c.neighbormolecules)) for c ∈ r.channels])
                                                                                                                                  for r ∈ rep]
                                                                                                                                  for rep ∈ sample]
                                                                                                                                  for sample ∈ experimentresults]
   
save(joinpath(outputdir, "resultsmols.jld2"), "experimentresults", experimentresults1)