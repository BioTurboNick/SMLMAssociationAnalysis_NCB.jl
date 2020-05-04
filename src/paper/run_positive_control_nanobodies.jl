# Recreates the analysis from the original data files

rootpath = "C:/Users/nicho/Dropbox (Partners HealthCare)/Data Analysis"
projectdirname = "MEG3 Project"
experimentdirnames = ["7 - U2OS FKBP12 mTOR STORM Nanobodies"]

datadirname = "Data"

samplenames = ["E", "G"]

nreplicates = 1
nsamples = 2
ncells = 10

outputdir = joinpath(rootpath, "SMLMAssociationAnalysis_NCB.jl", "original", "control nanobodies", "output")
mkpath(outputdir)
outputdatapath = joinpath(outputdir, "results.jld2")

using Distributed
currentworkers = addprocs(exeflags = "--project")
@everywhere using SMLMAssociationAnalysis_NCB
using Printf
using FileIO
using LocalizationMicroscopy
using Statistics


experimentresults = Vector{Vector{Vector{Result}}}[]
for experimentdirname ∈ experimentdirnames
    println("Starting experiment $experimentdirname.")
    experimentpath = joinpath(rootpath, projectdirname, experimentdirname, datadirname)
    replicateresults = Vector{Vector{Result}}[]
    for i ∈ 1:nreplicates
        sampleresults = Vector{Result}[]
        println("    Starting replicate $i.")
        replicatepath = joinpath(experimentpath, "Replicate $i")

        samplenames = i == 1 ? samplenames1 : samplenames2

        for samplename ∈ samplenames
            results = Result[]
            println("        Starting sample $samplename.")
            for j ∈ 1:ncells
                println("            Starting cell $j.")

                ch1_startframe = 1
                ch2_startframe = 11001

                cellpath = joinpath(replicatepath, "$samplename $(Printf.@sprintf("%03i", j)).bin.txt")
                localizations = loadlocalizations(cellpath, LocalizationMicroscopy.nikonelementstext)

                ch1_name = "647"
                ch2_name = "488"

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
                    10000
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

    save(outputdatapath, "replicateresults", experimentresults)
end

rmprocs(currentworkers)
