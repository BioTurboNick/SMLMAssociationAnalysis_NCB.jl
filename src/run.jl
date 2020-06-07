############
# PARAMETERS

# name of the project to save with the data
projectname = "My Project"

# the file paths of each file to analyze
# you may wish to write a loop to populate this
# The / path separator works on both Windows and Linux
cellpaths = ["paths/to/your/data/here1.bin.txt",
             "paths/to/your/data/here2.bin.txt"]

# jld2 is a julia-native file format. You can save to other formats by changing the extension and importing an appropriate package.
outputpath = "path/to/your/output/results.jld2"

# name of each channel
ch1_name = "647"
ch2_name = "488"

# number of the first frame in the acquistion
ch1_startframe = 1
ch2_startframe = 11001

# number of frames in the acquisition
ch1_frames = 11000
ch2_frames = 11000

# number of frames to remove from the start of acquisition
ch1_starttrim = 100
ch2_starttrim = 100

# number of frames to remove from the end of acquisition
ch1_endtrim = 10
ch2_endtrim = 10

# the maximum expected movement in mesured fluorophore position (nm)
maximum_displacement = 34.2

# the number of frames on either side of a localization to group within
t_off = 500

# the radius within which to merge localizations (nm)
merge_radius = 200

# Monte Carlo iterations to run
mc_iterations = 10000

# END PARAMETERS
################

using Distributed
currentworkers = addprocs(exeflags="--project")
@everywhere using SMLMAssociationAnalysis_NCB
using Statistics
using Printf
using LocalizationMicroscopy
using FileIO

results = Result[]

for i ∈ 1:length(cellpaths)
    cellpath = cellpaths[i]

    println("Processing $cellpath")
    localizations = loadlocalizations(cellpath, nikonelementstext)

    ch1_molecules, ch1_localizations = getmolecules(localizations, ch1_name, ch1_startframe, ch1_frames, ch1_starttrim,
                                                    ch1_endtrim, maximum_displacement, t_off, merge_radius)
    ch2_molecules, ch2_localizations = getmolecules(localizations, ch2_name, ch2_startframe, ch2_frames, ch2_starttrim,
                                                    ch2_endtrim, maximum_displacement, t_off, merge_radius)
    ch1_neighbors, ch2_neighbors, distances = exclusivenearestneighbors(ch1_molecules, ch2_molecules)

    percentileranks = montecarloaffinity(ch1_molecules, ch2_molecules, ch1_neighbors, ch2_neighbors, distances, 800, mc_iterations)
    mediandistance = length(distances) > 0 ? median(distances) : NaN

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

    ch1_data = ChannelData(ch1_name, ch1_molecules, ch1_neighbors)
    ch2_data = ChannelData(ch2_name, ch2_molecules, ch2_neighbors)
    result = Result(projectname, "", 1, cellpath, i,
                    [ch1_data, ch2_data], distances, mediandistance, percentileranks,
                    positivecontrol_distances, negativecontrol_distances,
                    positivecontrol_percentileranks, negativecontrol_percentileranks)
    push!(results, result)
end

mkpath(outputpath |> dirname)
save(outputpath, "results", results)

rmprocs(currentworkers)
