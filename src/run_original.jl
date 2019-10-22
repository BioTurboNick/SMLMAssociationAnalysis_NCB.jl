# Recreates the analysis from the original data files

datapath = raw"C:\Users\nicho\Dropbox (Partners HealthCare)\Data Analysis"
projectdirname = "MEG3 Project"
experimentdirnames = ["2 - U2OS p53 MDM2 STORM", "3 - U2OS p53 MEG3 STORM"]

datadirname = "Data"
outputdirname = "Output"

samplenames = ["A", "B", "C", "D"]

nreplicates = 3
nsamples = 4
ncells = 10

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
    experimentpath = joinpath(datapath, projectdirname, experimentdirname, datadirname)
    experimentoutputpath = joinpath(datapath, projectdirname, experimentdirname, outputdirname)
    replicateresults = Vector{Vector{Result}}[]
    for i ∈ 1:nreplicates
        sampleresults = Vector{Result}[]
        println("    Starting replicate $i.")
        replicatepath = joinpath(experimentpath, "Replicate $i")
        for samplename ∈ samplenames
            results = Result[]
            println("        Starting sample $samplename.")
            for j ∈ 1:ncells
                println("            Starting cell $j.")
                cellpath = joinpath(replicatepath, "$samplename $(Printf.@sprintf("%03i", j)).bin.txt")
                localizations = LocalizationMicroscopy.load(cellpath, LocalizationMicroscopy.nikonelementstext)
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
                    ch2_startframe = 11001 # try weeding out molecules with only 1 loc? ##############
                end
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

                percentileranks = montecarloaffinity(
                    ch1_molecules,
                    ch2_molecules,
                    ch1_neighbors,
                    ch2_neighbors,
                    distances,
                    200,
                    4,
                )

                if length(distances) == 0
                    mediandistance = NaN
                else
                    mediandistance = median(distances)
                end
                println("                Done: $(length(distances)) neighbors from $(length(ch1_molecules)) and $(length(ch2_molecules)) molecules, $(length(ch1_localizations)) and $(length(ch2_localizations)) localizations; median distance $mediandistance")
                ch1_data = ChannelData(ch1_name, ch1_molecules, ch1_neighbors)
                ch2_data = ChannelData(ch2_name, ch2_molecules, ch2_neighbors)
                result = Result(
                    projectdirname,
                    experimentdirname,
                    i,
                    samplename,
                    j,
                    [ch1_data, ch2_data],
                    distances,
                    mediandistance,
                    percentileranks,
                )
                push!(results, result)
            end
            push!(sampleresults, results)
        end
        push!(replicateresults, sampleresults)
    end
    push!(experimentresults, replicateresults)
end

outputpath = joinpath(experimentoutputpath, "results.jld2")
save(outputpath, "experimentresults", experimentresults)

rmprocs(currentworkers)
