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
#         - m = n = { 100, 200, 500, 1000, 2000 }
#         - d = { 10, 20, 50, 100, 200 } nm

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
moleculecounts = [100, 200, 500, 1000, 2000]
boundradii = [10, 20, 50, 100, 200]

cellradius = 2821 # 25 square micron circle

rootpath = "C:/Users/nicho/Dropbox (Partners HealthCare)/STORM MATLAB/STORM Single Molecule Clustering/MonteCarloAffinity/Simulated data"

using Distributed
using Tau
using Statistics
using LocalizationMicroscopy
using FileIO
currentworkers = addprocs(exeflags="--project")
@everywhere using SMLMAssociationAnalysis_NCB

results = Result[]

function generateunboundmolecules(molecules::Vector{Molecule}, moleculecount, cellradius)
    while length(molecules) < moleculecount
        unboundcoordinates = randomcoordinates2d(moleculecount - length(molecules), cellradius)
        unboundmolecules = [Molecule(Localization(length(molecules) + i, "", unboundcoordinates[1,i], unboundcoordinates[2,i], unboundcoordinates[3,i], 0, 1, 1)) for i ∈ 1:size(unboundcoordinates, 2)]
        append!(molecules, unboundmolecules)
        molecules = merge_close_molecules(molecules, 200)
    end
    molecules
end

for i ∈ 1:10
    for moleculecount ∈ moleculecounts
        for fractionbound ∈ fractionsbound
            boundcount = Int(round(fractionbound * moleculecount))

            for boundradius ∈ boundradii
                fractionbound == 0.0 && boundradius != boundradii[1] && break

                println("Processing $moleculecount mols $fractionbound bound $boundradius radius $i")

                boundcoordinates1 = randomcoordinates2d(boundcount, cellradius)
                molecules1 = [Molecule(Localization(i, "", boundcoordinates1[1,i], boundcoordinates1[2,i], boundcoordinates1[3,i], 0, 1, 1)) for i ∈ 1:size(boundcoordinates1, 2)]

                boundcoordinates2 = boundcoordinates1 .+ randomcoordinates2d(boundcount, boundradius)
                molecules2 = [Molecule(Localization(i, "", boundcoordinates2[1,i], boundcoordinates2[2,i], boundcoordinates2[3,i], 0, 1, 1)) for i ∈ 1:size(boundcoordinates2, 2)]

                molecules1 = generateunboundmolecules(molecules1, moleculecount, cellradius)
                molecules2 = generateunboundmolecules(molecules2, moleculecount, cellradius)

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
    save(joinpath(rootpath, "simulationresults$i.jld2"), "results", results)
end

rmprocs(currentworkers)


# plot

using StatsPlots

results = load(joinpath(rootpath, "simulationresults1.jld2"))["results"]

## 100 molecules each
# 0.0
distanceprobabilityplot(results[1])

# 0.01
distanceprobabilityplot(results[2])
distanceprobabilityplot(results[3])
distanceprobabilityplot(results[4])

# 0.02
distanceprobabilityplot(results[5])
distanceprobabilityplot(results[6])
distanceprobabilityplot(results[7])

# 0.05
distanceprobabilityplot(results[8])
distanceprobabilityplot(results[9])
distanceprobabilityplot(results[10])

# 0.10
distanceprobabilityplot(results[11])
distanceprobabilityplot(results[12])
distanceprobabilityplot(results[13])

# 0.20
distanceprobabilityplot(results[14])
distanceprobabilityplot(results[15])
distanceprobabilityplot(results[16])

# 0.50
distanceprobabilityplot(results[17])
distanceprobabilityplot(results[18])
distanceprobabilityplot(results[19])

# 1.0
distanceprobabilityplot(results[20])
distanceprobabilityplot(results[21])
distanceprobabilityplot(results[22])

moleculesplot_sim(results[6])


## 200 molecules each
# 0.0
distanceprobabilityplot(results[23])

# 0.01
distanceprobabilityplot(results[24])
distanceprobabilityplot(results[25])
distanceprobabilityplot(results[26])

# 0.02
distanceprobabilityplot(results[27])
distanceprobabilityplot(results[28])
distanceprobabilityplot(results[29])

# 0.05
distanceprobabilityplot(results[30])
distanceprobabilityplot(results[31])
distanceprobabilityplot(results[32])

# 0.10
distanceprobabilityplot(results[33])
distanceprobabilityplot(results[34])
distanceprobabilityplot(results[35])

# 0.20
distanceprobabilityplot(results[36])
distanceprobabilityplot(results[37])
distanceprobabilityplot(results[38])

# 0.50
distanceprobabilityplot(results[39])
distanceprobabilityplot(results[40])
distanceprobabilityplot(results[41])

# 1.0
distanceprobabilityplot(results[42])
distanceprobabilityplot(results[43])
distanceprobabilityplot(results[44])

## 500 molecules each
# 0.0
distanceprobabilityplot(results[45])

# 0.01
distanceprobabilityplot(results[46])
distanceprobabilityplot(results[47])
distanceprobabilityplot(results[48])

# 0.02
distanceprobabilityplot(results[49])
distanceprobabilityplot(results[50])
distanceprobabilityplot(results[51])

# 0.05
distanceprobabilityplot(results[52])
distanceprobabilityplot(results[53])
distanceprobabilityplot(results[54])

# 0.10
distanceprobabilityplot(results[55])
distanceprobabilityplot(results[56])
distanceprobabilityplot(results[57])

# 0.20
distanceprobabilityplot(results[58])
distanceprobabilityplot(results[59])
distanceprobabilityplot(results[60])

# 0.50
distanceprobabilityplot(results[61])
distanceprobabilityplot(results[62])
distanceprobabilityplot(results[63])

# 1.0
distanceprobabilityplot(results[64])
distanceprobabilityplot(results[65])
distanceprobabilityplot(results[66])
