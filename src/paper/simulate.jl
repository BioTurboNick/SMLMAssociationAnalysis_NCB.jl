# Simulate a cell:
#    1. Randomly generate positions of two sets of molecules of given amounts m, n
#         - select b = x% of min(m, n) binding pairs, for b pairs, m - b and n - b independent molcules
#         - randomly generate the m - b and n - b molecules
#         - generate the positions of bound pairs d distance apart, then randomly generate a 3D position and squash it
#           to 2D
#             - accounts for antibody distance too
#    2. Randomly generate blinks from each molecule (Poisson?)
#         - accounts for multiple labeling too
#    3. Run the simulated cell through the algorithm and obtain median distance and fraction bound
#    4. Repeat for varying combinations of m, n, x, d
#         - x = { 0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0 }
#         - m = { 100, 200, 500, 1000, 2000, 5000, 10000 },
#           n = { 100, 200, 500, 1000, 2000, 5000, 10000 }, m ≥ n
#         - d = { 20, 80, 140 } nm

# Basis for choices:
#
# In a particularly dense image, I have observed 2430 molecules and 500 molecules in 25 square microns, and in a less-
# dense area, 301 and 55 molecules, respectively.
#
# 20 nm  = 10 nm distance + 5 nm nanobodies
# 80 nm  = 70 nm distance + 5 nm nanobodies OR
#          10 nm distance + 35 nm antibody stack
# 140 nm = 70 nm distance + 35 nm antibody stack
#
# For context, microtubulues are 24 nm across, ribosomes 30 nm across
#

fractionsbound = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
moleculecounts = [100, 200, 500, 1000, 2000, 5000, 10000]
boundradii = [20, 80, 140]

cellradius = 2821 # 25 square micron circle

rootpath = "C:/Users/nicho/Dropbox (Partners HealthCare)/STORM MATLAB/STORM Single Molecule Clustering/MonteCarloAffinity/Simulated data"

using Distributed
using Tau
using Statistics
using FileIO
currentworkers = addprocs(exeflags="--project")
@everywhere using SMLMAssociationAnalysis_NCB

results = Result[]

for i ∈ 1:10
    for moleculecount ∈ moleculecounts
        for fractionbound ∈ fractionsbound
            boundcount = Int(round(fractionbound * moleculecount))
            unboundcount = moleculecount - boundcount

            unboundcoordinates1 = randomcoordinates2d(unboundcount, cellradius)
            unboundcoordinates2 = randomcoordinates2d(unboundcount, cellradius)
            boundcoordinates1 = randomcoordinates2d(boundcount, cellradius)

            for boundradius ∈ boundradii
                if fractionbound == 0.0 && boundradius != boundradii[1]
                    break
                end
                println("Processing $moleculecount mols $fractionbound bound $boundradius radius $i")

                boundcoordinates2 = boundcoordinates1 .+ randomcoordinates2d(boundcount, boundradius)

                coordinates1 = [unboundcoordinates1 boundcoordinates1]
                coordinates2 = [unboundcoordinates2 boundcoordinates2]

                molecules1 = [Molecule(Localization(i, "", coordinates1[1,i], coordinates1[2,i], coordinates1[3,i], 0, 1, 1)) for i ∈ 1:size(coordinates1, 2)]
                molecules2 = [Molecule(Localization(i, "", coordinates2[1,i], coordinates2[2,i], coordinates2[3,i], 0, 1, 1)) for i ∈ 1:size(coordinates2, 2)]

                neighbors1, neighbors2, distances = exclusivenearestneighbors(molecules1, molecules2)

                percentileranks = montecarloaffinity(molecules1, molecules2, neighbors1, neighbors2, distances, 200, 4, 10000)
                mediandistance = length(distances) > 0 ? median(distances) : NaN

                data1 = ChannelData("1", molecules1, neighbors1)
                data2 = ChannelData("2", molecules2, neighbors2)
                result = Result("", "", 1, "$moleculecount mols $fractionbound bound $boundradius radius", i,
                                [data1, data2], distances, mediandistance, percentileranks)
                push!(results, result)
            end
        end
    end
end

save(joinpath(rootpath, "simulationresults.jld2"), "results", results)

rmprocs(currentworkers)
