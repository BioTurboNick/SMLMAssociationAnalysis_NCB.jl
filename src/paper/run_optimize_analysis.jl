# Recreates the analysis from the original data files

using SMLMAssociationAnalysis_NCB
using FileIO
using StatsBase
using StatsPlots

### Load saved data

rootpath = "C:/Users/nicho/Dropbox (Partners HealthCare)/Data Analysis"

samplenames = ["A", "B", "C", "D"]

nreplicates = 3
nsamples = 4
ncells = 10

outputdir = joinpath(rootpath, "SMLMAssociationAnalysis_NCB.jl", "original", "output")
datapath = joinpath(outputdir, "results_optimize.jld2")

experimentresults = FileIO.load(datapath)["experimentresults"]

medianmeasurements = Array{Float64,4}(undef, ncells, nreplicates, 4, 2)
montecarlomeasurements = Array{Float64,4}(undef, ncells, nreplicates, 4, 2)
localization1counts = Array{Float64,4}(undef, ncells, nreplicates, 4, 2)
molecule1counts = Array{Float64,4}(undef, ncells, nreplicates, 4, 2)
localization2counts = Array{Float64,4}(undef, ncells, nreplicates, 4, 2)
molecule2counts = Array{Float64,4}(undef, ncells, nreplicates, 4, 2)

for k ∈ 1:2
    for i ∈ 1:nsamples
        samplemedianresults = Array{Float64,2}(undef, ncells, nreplicates)
        samplelessthan10results = Array{Float64,2}(undef, ncells, nreplicates)
        samplelocalization1counts = Array{Float64,2}(undef, ncells, nreplicates)
        samplemolecule1counts = Array{Float64,2}(undef, ncells, nreplicates)
        samplelocalization2counts = Array{Float64,2}(undef, ncells, nreplicates)
        samplemolecule2counts = Array{Float64,2}(undef, ncells, nreplicates)
        for j ∈ 1:nreplicates
            replicateresults = experimentresults[k][j][i]
            lessthanlimitreplicate = [(x.distances .< 200) .& (x.percentileranks .< 0.1) for x ∈ replicateresults]
            lessthan10 = count.(lessthanlimitreplicate) ./ length.(lessthanlimitreplicate)
            mediandistances = map(x -> x.mediandistance, replicateresults)
            samplemedianresults[:, j] = mediandistances
            samplelessthan10results[:, j] = lessthan10
            samplelocalization1counts[:, j] = [map(x -> x.group.localizations |> length, y.channels[1].molecules) |> sum for y ∈ replicateresults]
            samplelocalization2counts[:, j] = [map(x -> x.group.localizations |> length, y.channels[2].molecules) |> sum for y ∈ replicateresults]
            samplemolecule1counts[:, j] = map(x -> x.channels[1].molecules |> length, replicateresults)
            samplemolecule2counts[:, j] = map(x -> x.channels[2].molecules |> length, replicateresults)
        end
        medianmeasurements[:, :, i, k] = samplemedianresults
        montecarlomeasurements[:, :, i, k] = samplelessthan10results
        localization1counts[:, :, i, k] = samplelocalization1counts
        molecule1counts[:, :, i, k] = samplemolecule1counts
        localization2counts[:, :, i, k] = samplelocalization2counts
        molecule2counts[:, :, i, k] = samplemolecule2counts
    end
end


### Following optimization requires running run_optimize.jl instead of run_original.jl

# robustness (local distance)
average_percentileranks = [[[mean(x) for x ∈ experimentresults[1][k][j][i].percentileranks_by_localdensity]
                                     for i ∈ eachindex(experimentresults[1][k][j])]
                                     for j ∈ eachindex(experimentresults[1][k])]
                                     for k ∈ eachindex(experimentresults[1])

plot([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000], average_percentileranks)

savefig(joinpath(outputdir, "(partial) choice of local area effect on percentile rank.png"))

# temporal cutoff (needs special data)

ch1_counts = [[[[length(x) for x ∈ experimentresults[1][k][j][i].channels[1].molecules]
                           for i ∈ eachindex(experimentresults[1][k][j])]
                           for j ∈ eachindex(experimentresults[1][k])]
                           for k ∈ eachindex(experimentresults[1])]

ch1_frac = [[[[x ./ ch1_counts[k][j][i][1] for x ∈ ch1_counts[k][j][i]]
                                       for i ∈ eachindex(experimentresults[1][k][j])]
                                       for j ∈ eachindex(experimentresults[1][k])]
                                       for k ∈ eachindex(experimentresults[1])]

plot([0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000], ch1_frac)

savefig(joinpath(outputdir, "(partial) t_off selection molecules 1.png"))

ch2_counts = [[[[length(x) for x ∈ experimentresults[1][k][j][i].channels[2].molecules]
                           for i ∈ eachindex(experimentresults[1][k][j])]
                           for j ∈ eachindex(experimentresults[1][k])]
                           for k ∈ eachindex(experimentresults[1])]

ch2_frac = [[[[x ./ ch2_counts[k][j][i][1] for x ∈ ch2_counts[k][j][i]]
                            for i ∈ eachindex(experimentresults[1][k][j])]
                            for j ∈ eachindex(experimentresults[1][k])]
                            for k ∈ eachindex(experimentresults[1])]

plot([0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000], ch2_frac)

savefig(joinpath(outputdir, "(partial) t_off selection molecules 2.png"))
