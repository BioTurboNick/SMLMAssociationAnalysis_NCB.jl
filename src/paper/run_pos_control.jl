# Recreates the analysis from the original data files

datapath = "dataset"
projectname = "MEG3"
experimentdirnames = ["FKBP12-mTOR"]

samplenames1 = ["E", "G"]
samplenames23 = ["E", "F"]
samplenames4 = ["C", "F"]
samplenames5 = ["A", "B"]

nreplicates = 5
ncells = 10

mc_iterations = 10000

outputdir = "output"
mkpath(outputdir)
outputdatapath = joinpath(outputdir, "resultspos.jld2")

using Distributed
currentworkers = addprocs(exeflags = "--project")
@everywhere using SMLMAssociationAnalysis_NCB
using Printf
using LocalizationMicroscopy
using Statistics
using FileIO

experimentresults = Vector{Vector{Vector{Result}}}[]
for experimentdirname ∈ experimentdirnames
    println("Starting experiment $experimentdirname.")
    experimentpath = joinpath(datapath, experimentdirname)
    replicateresults = Vector{Vector{Result}}[]
    for i ∈ 1:nreplicates
        sampleresults = Vector{Result}[]
        println("    Starting replicate $i.")
        replicatepath = joinpath(experimentpath, "Replicate $i")
        samplenames = i == 5 ? samplenames5 :
                      i == 4 ? samplenames4 :
                      i == 1 ? samplenames1 :
                               samplenames23
        for samplename ∈ samplenames
            results = Result[]
            println("        Starting sample $samplename.")
            for j ∈ 1:ncells
                jcell = j
                if i == 1 && samplename == "E"
                    jcell += 10
                end

                println("            Starting cell $jcell.")
                cellpath = joinpath(replicatepath, "$samplename $(Printf.@sprintf("%03i", jcell)).bin.txt")
                localizations = loadlocalizations(cellpath, LocalizationMicroscopy.nikonelementstext)
                ch1_name = i == 2 ? "488" : "647"
                ch2_name = i == 2 ? "647" : "488"  # the channels were reversed for replicate 2. Unsure if this caused other issues.
                
                ch1_startframe = 1
                ch2_startframe = 11001
                
                ch1_task = remotecall(getmolecules, currentworkers[1], localizations,
                    ch1_name,
                    ch1_startframe,
                    11000,
                    100,
                    10,
                    34.2,
                    500,
                    200,
                )
                ch2_task = remotecall(getmolecules, currentworkers[2], localizations,
                    ch2_name,
                    ch2_startframe,
                    11000,
                    100,
                    10,
                    34.2,
                    500,
                    200,
                )

                if i == 2
                    ch1_task, ch2_task = ch2_task, ch1_task
                end

                ch1_molecules, ch1_localizations = fetch(ch1_task)
                ch2_molecules, ch2_localizations = fetch(ch2_task)

                ch1_neighbors, ch2_neighbors, distances = exclusivenearestneighbors(ch1_molecules, ch2_molecules)

                percentileranks = montecarloaffinity(
                    ch1_molecules,
                    ch2_molecules,
                    ch1_neighbors,
                    ch2_neighbors,
                    distances,
                    800,
                    mc_iterations
                )

                positivecontrol_percentileranks, positivecontrol_distances = simulate100(ch1_molecules, ch2_molecules, ch1_neighbors, ch2_neighbors, 80, 800, mc_iterations)
                # I found that the variance in running the negative control drops smoothly as a function of neighbor counts.
                # This heuristic cuts down on computational time without sacrificing quality of the correction.
                negative_iterations = length(distances) < 250 ? 30 :
                                      length(distances) < 500 ? 20 :
                                      length(distances) < 1000 ? 10 :
                                      length(distances) < 2000 ? 5 :
                                      1
                negativecontrol = [simulate0(ch1_molecules, ch2_molecules, 800, 10000) for i in 1:negative_iterations]
                negativecontrol_percentileranks = [first(i) for i in negativecontrol]
                negativecontrol_distances = [last(i) for i in negativecontrol]

                mediandistance = length(distances) == 0 ? NaN : median(distances)

                println("                Done: $(length(distances)) neighbors from $(length(ch1_molecules)) and $(length(ch2_molecules)) molecules, $(length(ch1_localizations)) and $(length(ch2_localizations)) localizations; median distance $mediandistance")
                ch1_data = ChannelData(ch1_name, ch1_molecules, ch1_neighbors)
                ch2_data = ChannelData(ch2_name, ch2_molecules, ch2_neighbors)
                result = Result(projectname, experimentdirname, 1, cellpath, i,
                                [ch1_data, ch2_data], distances, mediandistance, percentileranks,
                                positivecontrol_distances, negativecontrol_distances,
                                positivecontrol_percentileranks, negativecontrol_percentileranks)

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
