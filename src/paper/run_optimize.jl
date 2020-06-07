# Recreates the analysis from the original data files with different parameters

datapath = "dataset"
experimentdirnames = ["Mdm2-p53", "MEG3-p53"]

samplenames = ["A", "B", "C", "D"]

nreplicates = 3
nsamples = 4
ncells = 10

outputdir = "output/optimize"
mkpath(outputdir)
outputdatapath = joinpath(outputdir, "results_optimize.jld2")

using Distributed
currentworkers = addprocs(exeflags = "--project")
@everywhere using SMLMAssociationAnalysis_NCB
using Printf
using LocalizationMicroscopy
using Statistics
using FileIO

experimentresults = Vector{Vector{Vector{ResultOptimizing}}}[]
for experimentdirname ∈ experimentdirnames
    println("Starting experiment $experimentdirname.")
    experimentpath = joinpath(datapath, experimentdirname)
    replicateresults = Vector{Vector{ResultOptimizing}}[]
    for i ∈ 1:nreplicates
        sampleresults = Vector{ResultOptimizing}[]
        println("    Starting replicate $i.")
        replicatepath = joinpath(experimentpath, "Replicate $i")
        for samplename ∈ samplenames
            results = ResultOptimizing[]
            println("        Starting sample $samplename.")
            for j ∈ 1:ncells
                println("            Starting cell $j.")
                cellpath = joinpath(replicatepath, "$samplename $(Printf.@sprintf("%03i", j)).bin.txt")
                localizations = loadlocalizations(cellpath, LocalizationMicroscopy.nikonelementstext)
                # account for variances in data collection
                if experimentdirname == experimentdirnames[2] && samplename ∈ samplenames[3:4]
                    ch1_name = "561"
                else
                    ch1_name = "647"
                end
                ch2_name = "488"
                if experimentdirname == experimentdirnames[2] &&
                   ((i == 3 && samplename ∈ samplenames[1:3]) || (i == 2 && samplename == samplenames[2]))
                    ch1_startframe = 1
                    ch2_startframe = 15001
                else
                    ch1_startframe = 1
                    ch2_startframe = 11001
                end
                ch1 = getmolecules.(
                    Ref(localizations),
                    Ref(ch1_name),
                    Ref(ch1_startframe),
                    11000,
                    100,
                    10,
                    34.2,
                    [0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000],
                    0
                )
                ch2 = getmolecules.(
                    Ref(localizations),
                    Ref(ch2_name),
                    Ref(ch2_startframe),
                    11000,
                    100,
                    10,
                    34.2,
                    [0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000],
                    0
                )

                ch1_molecules, ch1_localizations = getmolecules(
                    localizations,
                    ch1_name,
                    ch1_startframe,
                    11000,
                    100,
                    10,
                    34.2,
                    500,
                    200,
                )
                ch2_molecules, ch2_localizations = getmolecules(
                    localizations,
                    ch2_name,
                    ch2_startframe,
                    11000,
                    100,
                    10,
                    34.2,
                    500,
                    200,
                )

                ch1_neighbors, ch2_neighbors, distances = exclusivenearestneighbors(ch1_molecules, ch2_molecules)

                percentileranks_by_localdensity = montecarloaffinity.(
                    Ref(ch1_molecules),
                    Ref(ch2_molecules),
                    Ref(ch1_neighbors),
                    Ref(ch2_neighbors),
                    Ref(distances),
                    200:200:2000,
                    10000
                )

                if length(distances) == 0
                    mediandistance = NaN
                else
                    mediandistance = median(distances)
                end
                println("                Done:$(length(distances)) neighbors from $(length(ch1_molecules)) and $(length(ch2_molecules)) molecules, $(length(ch1_localizations)) and $(length(ch2_localizations)) localizations; median distance $mediandistance")
                ch1_data = ChannelDataOptimizing(ch1_name, first.(ch1), ch1_neighbors)
                ch2_data = ChannelDataOptimizing(ch2_name, first.(ch2), ch2_neighbors)
                result = ResultOptimizing(
                    "",
                    experimentdirname,
                    i,
                    samplename,
                    j,
                    [ch1_data, ch2_data],
                    distances,
                    mediandistance,
                    percentileranks_by_localdensity
                )
                push!(results, result)
            end
            push!(sampleresults, results)
        end
        push!(replicateresults, sampleresults)
    end
    push!(experimentresults, replicateresults)
end

save(outputdatapath, "experimentresults", experimentresults)

rmprocs(currentworkers)
