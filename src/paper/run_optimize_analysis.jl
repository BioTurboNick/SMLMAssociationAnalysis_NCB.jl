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


### Following optimization requires running run_optimize.jl instead of run_original.jl
using Plots
using StatsBase

# robustness (local distance)
average_percentileranks = [[[[[mean(x) for x ∈ experimentresults[l][k][j][i].percentileranks_by_localdensity]
                                       for i ∈ eachindex(experimentresults[l][k][j])]
                                       for j ∈ eachindex(experimentresults[l][k])]
                                       for k ∈ eachindex(experimentresults[l])]
                                       for l ∈ eachindex(experimentresults)]

average_percentileranks = [[average_percentileranks[l][k][j][i][m] for i in 1:10
                                                                   for j in 1:4
                                                                   for k in 1:3
                                                                   for l in 1:2]
                                                                   for m in 1:10]

mean_average_percentileranks = mean.(average_percentileranks1)
sem_average_percentileranks = sem.(average_percentileranks1)
ci95_average_percentileranks = 1.96 .* sem_average_percentileranks

xvals = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
import Plots.mm
plot(xvals, average_percentileranks, line=(0.25))
plot!(xvals, mean_average_percentileranks, line=(stroke(10), :black), yaxis=(0:1), xaxis=(0:400:2000), tickfontsize=72, legend=:none, size=(1024,2048), left_margin=15mm, top_margin=5mm, right_margin=3mm, bottom_margin=5mm, ywiden=false)
plot!(xvals, mean_average_percentileranks - ci95_average_percentileranks, line=(stroke(10), :gray))
plot!(xvals, mean_average_percentileranks + ci95_average_percentileranks, line=(stroke(10), :gray))

savefig(joinpath(outputdir, "choice of local area effect on percentile rank.png"))

# temporal cutoff (needs special data)

ch1_counts = [[[[[length(x) for x ∈ experimentresults[l][k][j][i].channels[1].molecules]
                            for i ∈ eachindex(experimentresults[l][k][j])]
                            for j ∈ eachindex(experimentresults[l][k])]
                            for k ∈ eachindex(experimentresults[l])]
                            for l ∈ eachindex(experimentresults)]

ch1_frac = [[[[[x ./ ch1_counts[l][k][j][i][1] for x ∈ ch1_counts[l][k][j][i]]
                                            for i ∈ eachindex(experimentresults[l][k][j])]
                                            for j ∈ eachindex(experimentresults[l][k])]
                                            for k ∈ eachindex(experimentresults[l])]
                                            for l ∈ eachindex(experimentresults)]

ch1_frac1 = [[ch1_frac[l][k][j][i][m] for i in 1:10
                                       for j in 1:4
                                       for k in 1:3
                                       for l in 1:2]
                                       for m in 1:10]

mean_ch1_frac = mean.(ch1_frac1)
sem_ch1_frac = sem.(ch1_frac1)
ci95_ch1_frac = 1.96 .* sem_ch1_frac

xvals = [0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]

plot(xvals, ch1_frac, line=(0.25))
plot!(xvals, mean_ch1_frac, line=(stroke(10), :black), yaxis=(0:1), xaxis=(0:1000:5000), tickfontsize=64, guidefontsize=72, legend=:none, size=(1024,2048), left_margin=15mm, top_margin=5mm, right_margin=5mm, bottom_margin=5mm, ywiden=false)
plot!(xvals, mean_ch1_frac + ci95_ch1_frac, line=(stroke(10), :gray))
plot!(xvals, mean_ch1_frac - ci95_ch1_frac, line=(stroke(10), :gray))

savefig(joinpath(outputdir, "choice of t_off selection on number of molecules 1.png"))

ch2_counts = [[[[[length(x) for x ∈ experimentresults[l][k][j][i].channels[2].molecules]
                            for i ∈ eachindex(experimentresults[l][k][j])]
                            for j ∈ eachindex(experimentresults[l][k])]
                            for k ∈ eachindex(experimentresults[l])]
                            for l ∈ eachindex(experimentresults)]

ch2_frac = [[[[[x ./ ch2_counts[l][k][j][i][1] for x ∈ ch2_counts[l][k][j][i]]
                                                for i ∈ eachindex(experimentresults[l][k][j])]
                                                for j ∈ eachindex(experimentresults[l][k])]
                                                for k ∈ eachindex(experimentresults[l])]
                                                for l ∈ eachindex(experimentresults)]

ch2_frac1 = [[ch2_frac[l][k][j][i][m] for i in 1:10
                                       for j in 1:4
                                       for k in 1:3
                                       for l in 1:2]
                                       for m in 1:10]

mean_ch2_frac = mean.(ch2_frac1)
sem_ch2_frac = sem.(ch2_frac1)
ci95_ch2_frac = 1.96 .* sem_ch2_frac


xvals = [0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]

plot(xvals, ch2_frac, line=(0.25))
plot!(xvals, mean_ch2_frac, line=(stroke(10), :black), yaxis=(0:1), xaxis=(0:1000:5000), tickfontsize=64, guidefontsize=72, legend=:none, size=(1024,2048), left_margin=15mm, top_margin=5mm, bottom_margin=3mm, right_margin=5mm, ywiden=false)
plot!(xvals, mean_ch2_frac + ci95_ch2_frac, line=(stroke(10), :gray))
plot!(xvals, mean_ch2_frac - ci95_ch2_frac, line=(stroke(10), :gray))

savefig(joinpath(outputdir, "choice of t_off selection on number of molecules 2.png"))
