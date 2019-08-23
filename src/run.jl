using SMLMAssociationAnalysis_NCB

# PARAMETERS

# the file paths of each file to analyze
cellpaths = [raw"paths\to\your\data\here1.bin.txt",
             raw"paths\to\your\data\here2.bin.txt"]

outputpath = raw"path\to\your\output\results.jld2" # jld2 is a julia-native file format. You can save to other formats by changing the extension and importing an appropriate package.

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

# END PARAMETERS

results = Result[]

for i ∈ 1:length(cellpaths)
    cellpath = cellpaths[i]

    println("Processing $cellpath")
    localizations = LocalizationMicroscopy.load(cellpath, LocalizationMicroscopy.nikonelementstext)

    ch1_molecules, ch1_localizations = StormAnalysis.getmolecules(localizations, ch1_name, ch1_startframe, ch1_frames,
                                                                  ch1_starttrim, ch1_endtrim, maximum_displacement,
                                                                  t_off, merge_radius)
    ch2_molecules, ch2_localizations = StormAnalysis.getmolecules(localizations, ch2_name, ch2_startframe, ch2_frames,
                                                                  ch2_starttrim, ch2_endtrim, maximum_displacement,
                                                                  t_off, merge_radius)
    ch1_neighbors, ch2_neighbors, distances = StormAssociation.exclusivenearestneighbors(ch1_molecules, ch2_molecules)

    percentileranks = montecarloaffinity(ch1_molecules, ch2_molecules, ch1_neighbors, ch2_neighbors, distances, 4)

    if length(distances) == 0
        mediandistance = NaN
    else
        mediandistance = median(distances)
    end

    ch1_data = StormData.ChannelData(ch1_name, ch1_molecules, ch1_neighbors)
    ch2_data = StormData.ChannelData(ch2_name, ch2_molecules, ch2_neighbors)
    result = StormData.Result("", "", 1, cellpath, i,
                              [ch1_data, ch2_data], distances, mediandistance, percentileranks)
    push!(results, result)
end

save(outputpath, "results", results)
