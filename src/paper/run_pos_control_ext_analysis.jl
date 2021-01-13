# Recreates the analysis from the original data files

using SMLMAssociationAnalysis_NCB
using FileIO
using StatsBase
using Statistics
using StatsPlots
using SimpleANOVA

### Load saved data

samplenames = ["A", "B", "C", "D", "E", "F", "G", "H"]

nslides = 2
nsamples = 8
ncells = 10

outputdir = "output"
datapath = joinpath(outputdir, "results_ext.jld2")

experimentresults = load(datapath)["experimentresults"]

medianmeasurements = Array{Float64,3}(undef, ncells, nslides, 8)
montecarlomeasurements = Array{Float64,3}(undef, ncells, nslides, 8)
positivecontrolmontecarlomeasurements = Array{Float64,3}(undef, ncells, nslides, 8)
negativecontrolmontecarlomeasurements = Array{Float64,3}(undef, ncells, nslides, 8)

for i ∈ 1:nsamples
    samplemedianresults = Array{Float64,2}(undef, ncells, nslides)
    samplelessthan10results = Array{Float64,2}(undef, ncells, nslides)
    samplepositivecontrollessthan10results = Array{Float64,2}(undef, ncells, nslides)
    samplenegativecontrollessthan10results = Array{Float64,2}(undef, ncells, nslides)

    for j ∈ 1:nslides
        replicateresults = experimentresults[1][j][i]
        lessthanlimitreplicate = [(x.distances .< 200) .& (x.percentileranks .< 0.1) for x ∈ replicateresults]
        positivecontrollessthanlimitreplicate = [(x.positivecontrol_distances .< 200) .& (x.positivecontrol_percentileranks .< 0.1) for x ∈ replicateresults]
        negativecontrollessthanlimitreplicate = median.([count.([(x.negativecontrol_distances[m] .< 200) .& (x.negativecontrol_percentileranks[m] .< 0.1) for m ∈ eachindex(x.negativecontrol_distances)]) for x ∈ replicateresults])
        lessthan10 = count.(lessthanlimitreplicate) ./ length.(lessthanlimitreplicate)
        positivecontrollessthan10 = count.(positivecontrollessthanlimitreplicate) ./ length.(lessthanlimitreplicate)
        negativecontrollessthan10 = negativecontrollessthanlimitreplicate ./ length.(lessthanlimitreplicate)
        mediandistances = map(x -> x.mediandistance, replicateresults)
        samplemedianresults[:, j] = mediandistances
        samplelessthan10results[:, j] = lessthan10
        samplepositivecontrollessthan10results[:, j] = positivecontrollessthan10
        samplenegativecontrollessthan10results[:, j] = negativecontrollessthan10
    end
    medianmeasurements[:, :, i] = samplemedianresults
    montecarlomeasurements[:, :, i] = samplelessthan10results
    positivecontrolmontecarlomeasurements[:, :, i] = samplepositivecontrollessthan10results
    negativecontrolmontecarlomeasurements[:, :, i] = samplenegativecontrollessthan10results
end

normalizedmontecarlomeasurements = (montecarlomeasurements .- negativecontrolmontecarlomeasurements) ./ (positivecontrolmontecarlomeasurements .- negativecontrolmontecarlomeasurements)

#### Analysis
show("Automatically executing this section isn't ideal, as results may not all print, although graphics will be saved. Intended to go through one-by-one as part of the analysis.")

# plot FKBP12-RAPTOR

plot(permutedims(dropdims(mean(normalizedmontecarlomeasurements[:,:,1:2], dims=1), dims=1)), xticks=(1:2, ["NT", "+Rap"]), labels=["1h" "24h"], yerror=permutedims(dropdims(std(normalizedmontecarlomeasurements[:,:,1:2], dims=1), dims=1)))
savefig(joinpath(outputdir, "FKBP12-RAPTOR ext.png"))

# plot FKBP12-mTOR

plot([0, 2, 5, 10, 20, 40], permutedims(dropdims(mean(normalizedmontecarlomeasurements[:,:,3:8], dims=1), dims=1)), labels=["1h" "24h"], yerror=permutedims(dropdims(std(normalizedmontecarlomeasurements[:,:,3:8], dims=1), dims=1)))
savefig(joinpath(outputdir, "FKBP12-mTOR ext.png"))

# combined with previous positive control data
datapath = joinpath(outputdir, "results.jld2")
nsamples = 4
nreplicates = 3

experimentresults1 = load(datapath)["experimentresults"]

montecarlomeasurements1 = Array{Float64,3}(undef, ncells, nreplicates, 4)
positivecontrolmontecarlomeasurements1 = Array{Float64,3}(undef, ncells, nreplicates, 4)
negativecontrolmontecarlomeasurements1 = Array{Float64,3}(undef, ncells, nreplicates, 4)

for i ∈ 1:nsamples
    samplemedianresults = Array{Float64,2}(undef, ncells, nreplicates)
    samplelessthan10results = Array{Float64,2}(undef, ncells, nreplicates)
    samplepositivecontrollessthan10results = Array{Float64,2}(undef, ncells, nreplicates)
    samplenegativecontrollessthan10results = Array{Float64,2}(undef, ncells, nreplicates)

    for j ∈ 1:nreplicates
        replicateresults = experimentresults1[1][j][i]
        lessthanlimitreplicate = [(x.distances .< 200) .& (x.percentileranks .< 0.1) for x ∈ replicateresults]
        positivecontrollessthanlimitreplicate = [(x.positivecontrol_distances .< 200) .& (x.positivecontrol_percentileranks .< 0.1) for x ∈ replicateresults]
        negativecontrollessthanlimitreplicate = median.([count.([(x.negativecontrol_distances[m] .< 200) .& (x.negativecontrol_percentileranks[m] .< 0.1) for m ∈ eachindex(x.negativecontrol_distances)]) for x ∈ replicateresults])
        lessthan10 = count.(lessthanlimitreplicate) ./ length.(lessthanlimitreplicate)
        positivecontrollessthan10 = count.(positivecontrollessthanlimitreplicate) ./ length.(lessthanlimitreplicate)
        negativecontrollessthan10 = negativecontrollessthanlimitreplicate ./ length.(lessthanlimitreplicate)
        samplelessthan10results[:, j] = lessthan10
        samplepositivecontrollessthan10results[:, j] = positivecontrollessthan10
        samplenegativecontrollessthan10results[:, j] = negativecontrollessthan10
    end
    montecarlomeasurements1[:, :, i] = samplelessthan10results
    positivecontrolmontecarlomeasurements1[:, :, i] = samplepositivecontrollessthan10results
    negativecontrolmontecarlomeasurements1[:, :, i] = samplenegativecontrollessthan10results
end

normalizedmontecarlomeasurements1 = (montecarlomeasurements1 .- negativecontrolmontecarlomeasurements1) ./ (positivecontrolmontecarlomeasurements1 .- negativecontrolmontecarlomeasurements1)

# reoganize first replicate because its results are in a different order (0, 20, 10, 40)
montecarlomeasurements1 = [cat(montecarlomeasurements1[:,1,1], montecarlomeasurements1[:,1,3], montecarlomeasurements1[:,1,2], montecarlomeasurements1[:,1,4], dims = 3) montecarlomeasurements1[:,2:3,:]]
normalizedmontecarlomeasurements1 = [cat(normalizedmontecarlomeasurements1[:,1,1], normalizedmontecarlomeasurements1[:,1,3], normalizedmontecarlomeasurements1[:,1,2], normalizedmontecarlomeasurements1[:,1,4], dims = 3) normalizedmontecarlomeasurements1[:,2:3,:]]

combinednormalizedmontecarlomeasurements = [normalizedmontecarlomeasurements1[:,:,1:2] normalizedmontecarlomeasurements[:,2:2,[3,6]]]
combinedmontecarlomeasurements = [montecarlomeasurements1[:,:,1:2] montecarlomeasurements[:,2:2,[3,6]]]


#### Analysis
show("Automatically executing this section isn't ideal, as results may not all print, although graphics will be saved. Intended to go through one-by-one as part of the analysis.")

### Monte Carlo

# Check for unusual cases
p1 = boxplot(combinednormalizedmontecarlomeasurements[:, :, 1], xaxis = ("Replicates", [1, 2, 3, 4]), yaxis = ("Fraction bound"))
p2 = boxplot(combinednormalizedmontecarlomeasurements[:, :, 2], xaxis = ("Replicates", [1, 2, 3, 4]), yaxis = ("Fraction bound"))

plot(p1, p2, layout = grid(1, 2), legend = :none, plot_title = "FKBP12-mTOR")
savefig(joinpath(outputdir, "FKBP12-mTOR boxplots montecarlo.png"))


#ZResid/ZPred plot and Levene's test
montecarloflat = [combinednormalizedmontecarlomeasurements[:, :, 1] combinednormalizedmontecarlomeasurements[:, :, 2]]
z = zscore(montecarloflat)
zpred = repeat(mean(z, dims = 1), 10)
zresid = z .- zpred
scatter(
    zresid,
    zpred,
    xaxis = ("Standardized Residual (ZResid)"),
    yaxis = ("Standardized Predicted Value (ZPred)"),
    legend = :none,
)
savefig(joinpath(outputdir, "FKBP12-mTOR zresid-zpred montecarlo.png"))

levene(montecarloflat)

# qqnorm, skewness, kurtosis
p = [qqnorm(zresid[:, i]) for i ∈ 1:8]
plot(p..., layout = grid(4, 2), legend = :none, plot_title = "FKBP12-mTOR qqnorm montecarlo")
savefig(joinpath(outputdir, "FKBP12-mTOR qqnorm montecarlo.png"))

[skewness(montecarloflat[:, i]) for i ∈ 1:8]
[kurtosis(montecarloflat[:, i]) for i ∈ 1:8]

montecarloresult = anova(
    combinednormalizedmontecarlomeasurements[:, Not(3), :],
    [nested],
    factornames = ["Replicate", "Rapamycin"],
)

import Plots.mm
fkbp12_mtor_montecarlo = [combinednormalizedmontecarlomeasurements[:,:,1] |> vec; combinednormalizedmontecarlomeasurements[:,:,2] |> vec]
fkbp12_mtor_montecarlo = [combinedmontecarlomeasurements[:,:,1] |> vec; combinedmontecarlomeasurements[:,:,2] |> vec]
groups = repeat([1,2], inner = 40)
boxplot(groups, fkbp12_mtor_montecarlo, outliers=false,
        guidefontsize=36,
        tickfontsize=36,
        legend=:none,
        left_margin=25mm,
        top_margin=5mm,
        bottom_margin=5mm,
        seriescolor=[:white :lightgray],
        line=(6, 0.75),
        size=(1024,2048),
        xaxis=("Rapamycin", (1:2, ["-Rap", "+Rap"])),
        yaxis=("Fraction bound", (-0.1,0.2)))
fkbp12_mtor_montecarlo_means = dropdims(mean(combinednormalizedmontecarlomeasurements, dims=1), dims=1)
fkbp12_mtor_montecarlo_medians = dropdims(median(combinednormalizedmontecarlomeasurements, dims=1), dims=1)
dotplot!([1,1,1,1,2,2,2,2], fkbp12_mtor_montecarlo_medians |> vec, mode = :none, label="", marker=(12, 0.75, :rect, repeat([:orange, :darkblue, :darkred, :darkgreen]), stroke(0)))
dotplot!(groups, fkbp12_mtor_montecarlo, mode = :density, label="", marker=(8, 0.5, repeat([:orange, :darkblue, :darkred, :darkgreen], inner=10), stroke(0)))

savefig(joinpath(outputdir, "fkbp12_mtor_montecarlo_boxplot.png"))
