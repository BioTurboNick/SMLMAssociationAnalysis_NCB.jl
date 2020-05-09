# Recreates the analysis from the original data files

datapath = "C:/Users/nicho/Dropbox (Partners HealthCare)/Data Analysis" #"dataset"
experimentdirnames = ["Mdm2-p53", "3 - U2OS p53 MEG3 STORM"]

samplenames = ["A", "B", "C", "D"]

nreplicates = 3
nsamples = 4
ncells = 10

outputdir = "output"
mkpath(outputdir)
outputdatapath = joinpath(outputdir, "results4.jld2")

using Distributed
currentworkers = addprocs(exeflags = "--project")
@everywhere using SMLMAssociationAnalysis_NCB
using Printf
using LocalizationMicroscopy
using Statistics
using FileIO

# TODO: need to adjust for running all 1-4 replicates of MEG3-p53 vs. 3 for p53-Mdm2

experimentresults = Vector{Vector{Vector{Result}}}[]
for experimentdirname ∈ [experimentdirnames[2]] #experimentdirnames
    println("Starting experiment $experimentdirname.")
    experimentpath = joinpath(datapath, projectdirname, experimentdirname, datadirname)
    #experimentoutputpath = joinpath(datapath, projectdirname, experimentdirname, outputdirname)
    replicateresults = Vector{Vector{Result}}[]
    for i ∈ 4:4 #1:nreplicates
        sampleresults = Vector{Result}[]
        println("    Starting replicate $i.")
        replicatepath = joinpath(experimentpath, "Replicate $i")
        for samplename ∈ samplenames
            results = Result[]
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
                    800,
                    10000
                )

                positivecontrol_percentileranks = simulate100(ch1_molecules, ch2_molecules, ch1_neighbors, ch2_neighbors, 80, mc_iterations)
                negativecontrol_percentileranks = simulate0(ch1_molecules, ch2_molecules, 800, mc_iterations)

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
                    positivecontrol_percentileranks,
                    negativecontrol_percentileranks
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
