# Recreates the analysis from the original data files

using SMLMAssociationAnalysis_NCB
using FileIO
using StatsBase
using Statistics
using StatsPlots

### Load saved data

samplenames = ["A", "B", "C", "D"]

nreplicates = 3
nsamples = 4
ncells = 10

rootpath = "dataset"
outputdir = "output"
datapath = joinpath(rootpath, outputdir, "results.jld2")

experimentresults = load(datapath)["experimentresults"]

medianmeasurements = Array{Float64,4}(undef, ncells, nreplicates, 4, 2)
montecarlomeasurements = Array{Float64,4}(undef, ncells, nreplicates, 4, 2)
positivecontrolmontecarlomeasurements = Array{Float64,4}(undef, ncells, nreplicates, 4, 2)
negativecontrolmontecarlomeasurements = Array{Float64,4}(undef, ncells, nreplicates, 4, 2)

for k ∈ 1:2
    for i ∈ 1:nsamples
        samplemedianresults = Array{Float64,2}(undef, ncells, nreplicates)
        samplelessthan10results = Array{Float64,2}(undef, ncells, nreplicates)
        samplepositivecontrollessthan10results = Array{Float64,2}(undef, ncells, nreplicates)
        samplenegativecontrollessthan10results = Array{Float64,2}(undef, ncells, nreplicates)

        for j ∈ 1:nreplicates
            replicateresults = experimentresults[k][j][i]
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
        medianmeasurements[:, :, i, k] = samplemedianresults
        montecarlomeasurements[:, :, i, k] = samplelessthan10results
        positivecontrolmontecarlomeasurements[:, :, i, k] = samplepositivecontrollessthan10results
        negativecontrolmontecarlomeasurements[:, :, i, k] = samplenegativecontrollessthan10results
    end
end

medianmeasurements = cat(
    cat(medianmeasurements[:, :, 1:2, 1], medianmeasurements[:, :, 3:4, 1], dims = 4),
    cat(medianmeasurements[:, :, 1:2, 2], medianmeasurements[:, :, 3:4, 2], dims = 4),
    dims = 5,
)
montecarlomeasurements = cat(
    cat(montecarlomeasurements[:, :, 1:2, 1], montecarlomeasurements[:, :, 3:4, 1], dims = 4),
    cat(montecarlomeasurements[:, :, 1:2, 2], montecarlomeasurements[:, :, 3:4, 2], dims = 4),
    dims = 5,
)
positivecontrolmontecarlomeasurements = cat(
    cat(positivecontrolmontecarlomeasurements[:, :, 1:2, 1], positivecontrolmontecarlomeasurements[:, :, 3:4, 1], dims = 4),
    cat(positivecontrolmontecarlomeasurements[:, :, 1:2, 2], positivecontrolmontecarlomeasurements[:, :, 3:4, 2], dims = 4),
    dims = 5,
)
negativecontrolmontecarlomeasurements = cat(
    cat(negativecontrolmontecarlomeasurements[:, :, 1:2, 1], negativecontrolmontecarlomeasurements[:, :, 3:4, 1], dims = 4),
    cat(negativecontrolmontecarlomeasurements[:, :, 1:2, 2], negativecontrolmontecarlomeasurements[:, :, 3:4, 2], dims = 4),
    dims = 5,
)

normalizedmontecarlomeasurements = (montecarlomeasurements .- negativecontrolmontecarlomeasurements) ./ (positivecontrolmontecarlomeasurements .- negativecontrolmontecarlomeasurements)


#### Analysis
show("Automatically executing this section isn't ideal, as results may not all print, although graphics will be saved. Intended to go through one-by-one as part of the analysis.")

### Median Exp 2

# Check for unusual cases
p1 = boxplot(medianmeasurements[:, :, 1, 1, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
p2 = boxplot(medianmeasurements[:, :, 2, 1, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
p3 = boxplot(medianmeasurements[:, :, 1, 2, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
p4 = boxplot(medianmeasurements[:, :, 2, 2, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))

plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "MDM2-p53")
savefig(joinpath(outputdir, "Mdm2-p53 boxplots.png"))


#ZResid/ZPred plot and Levene's test
mediansflat = [medianmeasurements[:, :, 1, 1, 1] medianmeasurements[:, :, 2, 1, 1] medianmeasurements[:, :, 1, 2, 1] medianmeasurements[
    :,
    :,
    2,
    2,
    1,
]]
z = zscore(mediansflat)
zpred = repeat(mean(z, dims = 1), 10)
zresid = z .- zpred
scatter(
    zresid,
    zpred,
    xaxis = ("Standardized Residual (ZResid)"),
    yaxis = ("Standardized Predicted Value (ZPred)"),
    legend = :none,
)
savefig(joinpath(outputdir, "Mdm2-p53 zresid-zpred.png"))

levene(mediansflat)

# qqnorm, skewness, kurtosis
p = [qqnorm(zresid[:, i]) for i ∈ 1:12]
plot(p..., layout = grid(4, 3), legend = :none, plot_title = "MDM2-p53 qqnorm")
savefig(joinpath(outputdir, "Mdm2-p53 qqnorm.png"))

[skewness(mediansflat[:, i]) for i ∈ 1:12]
[kurtosis(mediansflat[:, i]) for i ∈ 1:12]


# anova

medianresult = anova(
    medianmeasurements[:, :, :, :, 1],
    [nested],
    factornames = ["Replicate", "Doxycycline", "Nutlin-3a"],
)

plot(medianresult)
savefig(joinpath(outputdir, "Mdm2-p53 interactionplot.png"))




### Monte Carlo Exp 2

# Check for unusual cases
p1 = boxplot(montecarlomeasurements[:, :, 1, 1, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p2 = boxplot(montecarlomeasurements[:, :, 2, 1, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p3 = boxplot(montecarlomeasurements[:, :, 1, 2, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p4 = boxplot(montecarlomeasurements[:, :, 2, 2, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))

plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "Mdm2-p53")
savefig(joinpath(outputdir, "Mdm2-p53 boxplots montecarlo.png"))


#ZResid/ZPred plot and Levene's test
montecarloflat = [montecarlomeasurements[:, :, 1, 1, 1] montecarlomeasurements[:, :, 2, 1, 1] montecarlomeasurements[:, :, 1, 2, 1] montecarlomeasurements[:, :, 2, 2, 1]]
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
savefig(joinpath(outputdir, "Mdm2-p53 zresid-zpred montecarlo.png"))

levene(montecarloflat)

# qqnorm, skewness, kurtosis
p = [qqnorm(zresid[:, i]) for i ∈ 1:12]
plot(p..., layout = grid(4, 3), legend = :none, plot_title = "MDM2-p53 qqnorm montecarlo")
savefig(joinpath(outputdir, "Mdm2-p53 qqnorm montecarlo.png"))

[skewness(montecarloflat[:, i]) for i ∈ 1:12]
[kurtosis(montecarloflat[:, i]) for i ∈ 1:12]


# anova

montecarloresult = anova(
    montecarlomeasurements[:, :, :, :, 1],
    [nested],
    factornames = ["Replicate", "Doxycycline", "Nutlin-3a"],
)

plot(montecarloresult)
savefig(joinpath(outputdir, "Mdm2-p53 interactionplot montecarlo.png"))

montecarloresultNutMinus = anova(
    montecarlomeasurements[:, :, :, 1, 2],
    [nested],
    factornames = ["Replicate", "Doxycycline"],
)
montecarloresultNutPlus = anova(
    montecarlomeasurements[:, :, :, 2, 2],
    [nested],
    factornames = ["Replicate", "Doxycycline"],
)

### Paper figure

nutdoxgroups = repeat([2, 1, 4, 3], inner=30)
p53_mdm2_median = [medianmeasurements[:,:,1,1,1] |> vec; medianmeasurements[:,:,2,1,1] |> vec; medianmeasurements[:,:,1,2,1] |> vec; medianmeasurements[:,:,2,2,1] |> vec]
boxplot(nutdoxgroups, p53_mdm2_median, outliers=false,
        label=["- Dox", "+ Dox"],
        guidefontsize=52,
        tickfontsize=48,
        legend=:none,
        top_margin=10mm,
        bottom_margin=-25mm,
        left_margin=-3mm,
        seriescolor=[:white :lightgray],
        line=(6, 1.0),
        gridopacity=0.3,
        xgrid=:none,
        size=(1024,2048),
        xaxis=("Condition", (1:4, ["-Dox -Nut", "+Dox -Nut", "-Dox +Nut", "+Dox +Nut"])),
        yaxis=("Median exclusive pairwise distance (nm)", (0,500)))
dotplot!(nutdoxgroups, p53_mdm2_median, mode = :density, label="", marker=(8, repeat([:orange, :darkblue, :darkred], inner=10), stroke(0)))
savefig( "p53_mdm2_median_boxplot.png")

p53_mdm2_montecarlo = [montecarlomeasurements[:,:,1,1,1] |> vec; montecarlomeasurements[:,:,2,1,1] |> vec; montecarlomeasurements[:,:,1,2,1] |> vec; montecarlomeasurements[:,:,2,2,1] |> vec]
groupedboxplot(nutdoxgroups, p53_mdm2_montecarlo, outliers=false,
        label=["- Dox", "+ Dox"],
        guidefontsize=12,
        tickfontsize=72,
        legend=:none,
        left_margin=20mm,
        top_margin=5mm,
        bottom_margin=5mm,
        seriescolor=[:white :lightgray],
        line=(6, 1.0),
        size=(1024,2048),
        xaxis=("Condition", (1:4, ["-Dox -Nut", "+Dox -Nut", "-Dox +Nut", "+Dox +Nut"])),
        yaxis=("Fraction bound", (0,0.15)))
groupeddotplot!(nutdoxgroups, p53_mdm2_montecarlo, mode = :density, label="", marker=(8, repeat([:orange, :darkblue, :darkred], inner=10), stroke(0)))
savefig(joinpath(outputdir, "p53_mdm2_montecarlo_boxplot.png"))

p53_mdm2_montecarlo = [normalizedmontecarlomeasurements[:,:,1,1,1] |> vec; normalizedmontecarlomeasurements[:,:,2,1,1] |> vec; normalizedmontecarlomeasurements[:,:,1,2,1] |> vec; normalizedmontecarlomeasurements[:,:,2,2,1] |> vec]
groupedboxplot(nutdoxgroups, p53_mdm2_montecarlo, outliers=false,
        label=["- Dox", "+ Dox"],
        guidefontsize=52,
        tickfontsize=48,
        legend=:none,
        top_margin=10mm,
        bottom_margin=-25mm,
        left_margin=5mm,
        seriescolor=[:white :lightgray],
        line=(6, 1.0),
        gridopacity=0.3,
        xgrid=:none,
        size=(1024,2048),
        xaxis=("Condition", (1:4, ["-Dox -Nut", "+Dox -Nut", "-Dox +Nut", "+Dox +Nut"])),
        yaxis=("Fraction bound", (-0.04,0.1), -0.04:0.02:0.1))
groupeddotplot!(nutdoxgroups, p53_mdm2_montecarlo, mode = :density, label="",
        marker=(8, repeat([:orange, :darkblue, :darkred], inner=10), stroke(0)))
savefig("p53_mdm2_normalizedmontecarlo_boxplot.png")

### Median Exp 3

# Check for unusual cases
p1 = boxplot(medianmeasurements[:, :, 1, 1, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
p2 = boxplot(medianmeasurements[:, :, 2, 1, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
p3 = boxplot(medianmeasurements[:, :, 1, 2, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
p4 = boxplot(medianmeasurements[:, :, 2, 2, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))

plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "MEG3-p53")
savefig(joinpath(outputdir, "MEG3-p53 boxplots.png"))


#ZResid/ZPred plot and Levene's test
mediansflat = [medianmeasurements[:, :, 1, 1, 2] medianmeasurements[:, :, 2, 1, 2] medianmeasurements[:, :, 1, 2, 2] medianmeasurements[:, :, 2, 2, 2]]
z = zscore(mediansflat)
zpred = repeat(mean(z, dims = 1), 10)
zresid = z .- zpred
scatter(
    zresid,
    zpred,
    xaxis = ("Standardized Residual (ZResid)"),
    yaxis = ("Standardized Predicted Value (ZPred)"),
    legend = :none,
)
savefig(joinpath(outputdir, "MEG3-p53 zresid-zpred.png"))

levene(mediansflat)

# qqnorm, skewness, kurtosis
p = [qqnorm(zresid[:, i]) for i ∈ 1:12]
plot(p..., layout = grid(4, 3), legend = :none, plot_title = "MEG3-p53 qqnorm")
savefig(joinpath(outputdir, "MEG3-p53 qqnorm.png"))

[skewness(mediansflat[:, i]) for i ∈ 1:12]
[kurtosis(mediansflat[:, i]) for i ∈ 1:12]


# anova

medianresult = anova(medianmeasurements[:, :, :, :, 2], [nested], factornames = ["Replicate", "Doxycycline", "RNA"])

plot(medianresult)
savefig(joinpath(outputdir, "MEG3-p53 interactionplot.png"))




### Monte Carlo Exp 3

# Check for unusual cases
p1 = boxplot(montecarlomeasurements[:, :, 1, 1, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p2 = boxplot(montecarlomeasurements[:, :, 2, 1, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p3 = boxplot(montecarlomeasurements[:, :, 1, 2, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p4 = boxplot(montecarlomeasurements[:, :, 2, 2, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))

plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "MEG3-p53")
savefig(joinpath(outputdir, "MEG3-p53 boxplots montecarlo.png"))


#ZResid/ZPred plot and Levene's test
montecarloflat = [montecarlomeasurements[:, :, 1, 1, 2] montecarlomeasurements[:, :, 2, 1, 2] montecarlomeasurements[:, :, 1, 2, 2] montecarlomeasurements[:, :, 2, 2, 2]]
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
savefig(joinpath(outputdir, "MEG3-p53 zresid-zpred montecarlo.png"))

levene(montecarloflat)

# qqnorm, skewness, kurtosis
p = [qqnorm(zresid[:, i]) for i ∈ 1:12]
plot(p..., layout = grid(4, 3), legend = :none, plot_title = "MEG3-p53 qqnorm montecarlo")
savefig(joinpath(outputdir, "MEG3-p53 qqnorm montecarlo.png"))

[skewness(montecarloflat[:, i]) for i ∈ 1:12]
[kurtosis(montecarloflat[:, i]) for i ∈ 1:12]


# anova

montecarloresult = anova(
    montecarlomeasurements[:, :, :, :, 2],
    [nested],
    factornames = ["Replicate", "Doxycycline", "RNA"],
)

plot(montecarloresult)
savefig(joinpath(outputdir, "MEG3-p53 interactionplot montecarlo.png"))

montecarloresultMEG3 = anova(
    montecarlomeasurements[:, :, :, 1, 2],
    [nested],
    factornames = ["Replicate", "Doxycycline"],
)
montecarloresultGAPDH = anova(
    montecarlomeasurements[:, :, :, 2, 2],
    [nested],
    factornames = ["Replicate", "Doxycycline"],
)

### Paper figure

import Plots.mm
rnagroups = repeat([1,2], inner = 60)
doxgroups = repeat([2, 1], inner = 30, outer = 2)
p53_meg3_median = [medianmeasurements[:,:,1,1,2] |> vec; medianmeasurements[:,:,2,1,2] |> vec; medianmeasurements[:,:,1,2,2] |> vec; medianmeasurements[:,:,2,2,2] |> vec]
groupedboxplot(rnagroups, p53_meg3_median, group = doxgroups, outliers=false,
        label=["- Dox", "+ Dox"],
        guidefontsize=52,
        tickfontsize=48,
        legend=:none,
        top_margin=10mm,
        bottom_margin=-25mm,
        seriescolor=[:white :lightgray],
        line=(6, 1.0),
        gridopacity=0.3,
        xgrid=:none,
        xaxis=("RNA", (1:2, ["MEG3", "GAPDH"])),
        size=(1024,2048),
        yaxis=("Median exclusive pairwise distance (nm)", (0,1200)))
groupeddotplot!(rnagroups, p53_meg3_median, group = doxgroups, mode = :density, label="", marker=(8, repeat([:orange, :darkblue, :darkred], inner=10), stroke(0)))
savefig("p53_meg3_median_boxplot.png")

p53_meg3_montecarlo = [montecarlomeasurements[:,:,1,1,2] |> vec; montecarlomeasurements[:,:,2,1,2] |> vec; montecarlomeasurements[:,:,1,2,2] |> vec; montecarlomeasurements[:,:,2,2,2] |> vec]
groupedboxplot(rnagroups, p53_meg3_montecarlo, group = doxgroups, outliers=false,
        label=["- Dox", "+ Dox"],
        guidefontsize=12,
        tickfontsize=72,
        legend=:none,
        left_margin=20mm,
        top_margin=5mm,
        bottom_margin=5mm,
        seriescolor=[:white :lightgray],
        line=(6, 1.0),
        size=(1024,2048),
        xaxis=("RNA", (1:2, ["MEG3", "GAPDH"])),
        yaxis=("Fraction bound", (0,0.2)))
groupeddotplot!(rnagroups, p53_meg3_montecarlo, group = doxgroups, mode = :density, label="", marker=(8, repeat([:orange, :darkblue, :darkred], inner=10), stroke(0)))
savefig(joinpath(outputdir, "p53_meg3_montecarlo_boxplot.png"))

p53_meg3_montecarlo = [normalizedmontecarlomeasurements[:,:,1,1,2] |> vec; normalizedmontecarlomeasurements[:,:,2,1,2] |> vec; normalizedmontecarlomeasurements[:,:,1,2,2] |> vec; normalizedmontecarlomeasurements[:,:,2,2,2] |> vec]
groupedboxplot(rnagroups, p53_meg3_montecarlo, group = doxgroups, outliers=false,
        label=["- Dox", "+ Dox"],
        guidefontsize=52,
        tickfontsize=48,
        legend=:none,
        top_margin=10mm,
        bottom_margin=-25mm,
        seriescolor=[:white :lightgray],
        line=(6, 1.0),
        gridopacity=0.3,
        xgrid=:none,
        size=(1024,2048),
        xaxis=("RNA", (1:2, ["MEG3", "GAPDH"])),
        yaxis=("Fraction bound", (-0.05,0.2), -0.05:0.05:0.2))
groupeddotplot!(rnagroups, p53_meg3_montecarlo, group = doxgroups, mode = :density, label="",
        marker=(8, repeat([:orange, :darkblue, :darkred], inner=10), stroke(0)))
savefig(joinpath(outputdir, "p53_meg3_normalizedmontecarlo_boxplot.png"))

# localization plots of example cells
# Exp 3 (p53-MEG3)
# A2
insetx, insety = [11000, 15096], [19000, 23096]
localizationsplot_forprint(experimentresults[2][1][1][2], insetbox = [insetx, insety])
savefig(joinpath(outputdir, "3 - A2 dSTORM points.png"))
localizationsinsetplot_forprint(experimentresults[2][1][1][2], insetx, insety)
savefig(joinpath(outputdir, "3 - A2 dSTORM points 10x.png"))
moleculesinsetplot_forprint(experimentresults[2][1][1][2], insetx, insety)
savefig(joinpath(outputdir, "3 - A2 dSTORM molecules 10x.png"))
insetplot_forprint(experimentresults[2][1][1][2], insetx, insety)
savefig(joinpath(outputdir, "3 - A2 dSTORM 10x.png"))

# B2
insetx, insety = [17000, 21096], [18000, 22096]
localizationsplot_forprint(experimentresults[2][1][2][2], insetbox = [insetx, insety])
savefig(joinpath(outputdir, "3 - B2 dSTORM points.png"))
localizationsinsetplot_forprint(experimentresults[2][1][2][2], insetx, insety)
savefig(joinpath(outputdir, "3 - B2 dSTORM points 10x.png"))
moleculesinsetplot_forprint(experimentresults[2][1][2][2], insetx, insety)
savefig(joinpath(outputdir, "3 - B2 dSTORM molecules 10x.png"))
insetplot_forprint(experimentresults[2][1][2][2], insetx, insety)
savefig(joinpath(outputdir, "3 - B2 dSTORM 10x.png"))

# C10
insetx, insety = [7500, 11596], [20000, 24096]
localizationsplot_forprint(experimentresults[2][1][3][10], insetbox = [insetx, insety])
savefig(joinpath(outputdir, "3 - C10 dSTORM points.png"))
localizationsinsetplot_forprint(experimentresults[2][1][3][10], insetx, insety)
savefig(joinpath(outputdir, "3 - C10 dSTORM points 10x.png"))
moleculesinsetplot_forprint(experimentresults[2][1][3][10], insetx, insety)
savefig(joinpath(outputdir, "3 - C10 dSTORM molecules 10x.png"))
insetplot_forprint(experimentresults[2][1][3][10], insetx, insety)
savefig(joinpath(outputdir, "3 - C10 dSTORM 10x.png"))


# D10
insetx, insety = [7500, 11596], [20000, 24096]
localizationsplot_forprint(experimentresults[2][1][4][10], insetbox = [insetx, insety])
savefig(joinpath(outputdir, "3 - D10 dSTORM points.png"))
localizationsinsetplot_forprint(experimentresults[2][1][4][10], insetx, insety)
savefig(joinpath(outputdir, "3 - D10 dSTORM points 10x.png"))
moleculesinsetplot_forprint(experimentresults[2][1][4][10], insetx, insety)
savefig(joinpath(outputdir, "3 - D10 dSTORM molecules 10x.png"))
insetplot_forprint(experimentresults[2][1][4][10], insetx, insety)
savefig(joinpath(outputdir, "3 - D10 dSTORM 10x.png"))

# Exp 2 (p53-MDM2)
# R3 A8
insetx, insety = [14000, 18096], [28000, 32096]
localizationsplot_forprint(experimentresults[1][3][1][8], insetbox = [insetx, insety])
savefig(joinpath(outputdir, "2 - A8 dSTORM points.png"))
localizationsinsetplot_forprint(experimentresults[1][3][1][8], insetx, insety)
savefig(joinpath(outputdir, "2 - A8 dSTORM points 10x.png"))
moleculesinsetplot_forprint(experimentresults[1][3][1][8], insetx, insety)
savefig(joinpath(outputdir, "2 - A8 dSTORM molecules 10x.png"))
insetplot_forprint(experimentresults[1][3][1][8], insetx, insety)
savefig(joinpath(outputdir, "2 - A8 dSTORM 10x.png"))

# R3 B5
insetx, insety = [14000, 18096], [20000, 24096]
localizationsplot_forprint(experimentresults[1][3][2][5], insetbox = [insetx, insety])
savefig(joinpath(outputdir, "2 - B5 dSTORM points.png"))
localizationsinsetplot_forprint(experimentresults[1][3][2][5], insetx, insety)
savefig(joinpath(outputdir, "2 - B5 dSTORM points 10x.png"))
moleculesinsetplot_forprint(experimentresults[1][3][2][5], insetx, insety)
savefig(joinpath(outputdir, "2 - B5 dSTORM molecules 10x.png"))
insetplot_forprint(experimentresults[1][3][2][5], insetx, insety)
savefig(joinpath(outputdir, "2 - B5 dSTORM 10x.png"))

# R2 C9
insetx, insety = [15000, 19096], [20000, 24096]
localizationsplot_forprint(experimentresults[1][2][3][9], insetbox = [insetx, insety])
savefig(joinpath(outputdir, "2 - C9 dSTORM points.png"))
localizationsinsetplot_forprint(experimentresults[1][2][3][9], insetx, insety)
savefig(joinpath(outputdir, "2 - C9 dSTORM points 10x.png"))
moleculesinsetplot_forprint(experimentresults[1][2][3][9], insetx, insety)
savefig(joinpath(outputdir, "2 - C9 dSTORM molecules 10x.png"))
insetplot_forprint(experimentresults[1][2][3][9], insetx, insety)
savefig(joinpath(outputdir, "2 - C9 dSTORM 10x.png"))

# R2 D4
insetx, insety = [21500, 25596], [12500, 16596]
localizationsplot_forprint(experimentresults[1][2][4][4], insetbox = [insetx, insety])
savefig(joinpath(outputdir, "2 - D4 dSTORM points.png"))
localizationsinsetplot_forprint(experimentresults[1][2][4][4], insetx, insety)
savefig(joinpath(outputdir, "2 - D4 dSTORM points 10x.png"))
moleculesinsetplot_forprint(experimentresults[1][2][4][4], insetx, insety)
savefig(joinpath(outputdir, "2 - D4 dSTORM molecules 10x.png"))
insetplot_forprint(experimentresults[1][2][4][4], insetx, insety)
savefig(joinpath(outputdir, "2 - D4 dSTORM 10x.png"))
