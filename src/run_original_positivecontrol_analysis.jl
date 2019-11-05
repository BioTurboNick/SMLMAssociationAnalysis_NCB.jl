# Recreates the analysis from the original data files

using SMLMAssociationAnalysis_NCB
using FileIO
using StatsBase
using StatsPlots

### Load saved data

rootpath = "C:/Users/nicho/Dropbox (Partners HealthCare)/Data Analysis"

samplenames = ["E", "F", "G", "H"]

nreplicates = 1
nsamples = 4
ncells = 10

outputdir = joinpath(rootpath, "SMLMAssociationAnalysis_NCB.jl", "original", "control", "output")
datapath = joinpath(outputdir, "results.jld2")

experimentresults = FileIO.load(datapath)["replicateresults"]

medianmeasurements = Array{Float64,3}(undef, ncells, nreplicates, 4)
montecarlomeasurements = Array{Float64,3}(undef, ncells, nreplicates, 4)

for k ∈ 1:1
    for i ∈ 1:nsamples
        samplemedianresults = Array{Float64,2}(undef, ncells, nreplicates)
        samplelessthan10results = Array{Float64,2}(undef, ncells, nreplicates)
        for j ∈ 1:nreplicates
            replicateresults = experimentresults[k][j][i]
            percentilerankslessthan10replicate = [x.percentileranks .< 0.1 for x ∈ replicateresults]
            lessthan10 = count.(percentilerankslessthan10replicate) ./ length.(percentilerankslessthan10replicate)
            mediandistances = map(x -> x.mediandistance, replicateresults)
            samplemedianresults[:, j] = mediandistances
            samplelessthan10results[:, j] = lessthan10
        end
        medianmeasurements[:, :, i] = samplemedianresults
        montecarlomeasurements[:, :, i] = samplelessthan10results
    end
end

medianmeasurements_temp = medianmeasurements[:, :, 2]
medianmeasurements[:, :, 2] = medianmeasurements[:, :, 3]
medianmeasurements[:, :, 3] = medianmeasurements_temp

montecarlomeasurements_temp = montecarlomeasurements[:, :, 2]
montecarlomeasurements[:, :, 2] = montecarlomeasurements[:, :, 3]
montecarlomeasurements[:, :, 3] = montecarlomeasurements_temp

#### Analysis
show("Automatically executing this section isn't ideal, as results may not all print, although graphics will be saved. Intended to go through one-by-one as part of the analysis.")

# Check for unusual cases
p1 = boxplot(medianmeasurements[:, :, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
p2 = boxplot(medianmeasurements[:, :, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
p3 = boxplot(medianmeasurements[:, :, 3], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
p4 = boxplot(medianmeasurements[:, :, 4], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))

plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "FKBP12-mTOR")


# Windsorize an extreme value
medianmeasurements[4, 1, 1] = medianmeasurements[5, 1, 1]

savefig(joinpath(outputdir, "FKBP12-mTOR boxplots-afterwindsorizing.png"))


#ZResid/ZPred plot and Levene's test
mediansflat = [medianmeasurements[:, :, 1] medianmeasurements[:, :, 2] medianmeasurements[:, :, 3] medianmeasurements[:, :, 4]]
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
savefig(joinpath(outputdir, "FKBP12-mTOR zresid-zpred.png"))

using SimpleANOVA
levene(mediansflat)

# qqnorm, skewness, kurtosis
using StatsPlots
p = [qqnorm(zresid[:, i]) for i ∈ 1:4]
plot(p..., layout = grid(4, 1), legend = :none, plot_title = "MEG3-p53 qqnorm")
savefig(joinpath(outputdir, "FKBP12-mTOR qqnorm.png"))

[skewness(mediansflat[:, i]) for i ∈ 1:4]
[kurtosis(mediansflat[:, i]) for i ∈ 1:4]


# anova

medianresult = anova(medianmeasurements[:, 1, :], factornames = ["Rapamycin"])
differencecontrasts(medianresult)

### Monte Carlo Exp 3

# Check for unusual cases
p1 = boxplot(montecarlomeasurements[:, :, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p2 = boxplot(montecarlomeasurements[:, :, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p3 = boxplot(montecarlomeasurements[:, :, 3], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p4 = boxplot(montecarlomeasurements[:, :, 4], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))

plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "FKBP12-mTOR")
savefig(joinpath(outputdir, "FKBP12-mTOR boxplots montecarlo.png"))

#ZResid/ZPred plot and Levene's test
montecarloflat = [montecarlomeasurements[:, :, 1] montecarlomeasurements[:, :, 2] montecarlomeasurements[:, :, 3] montecarlomeasurements[:, :, 4]]
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

using SimpleANOVA
levene(montecarloflat)

# qqnorm, skewness, kurtosis
using StatsPlots
p = [qqnorm(zresid[:, i]) for i ∈ 1:4]
plot(p..., layout = grid(4, 1), legend = :none, plot_title = "FKBP12-mTOR qqnorm montecarlo")
savefig(joinpath(outputdir, "FKBP12-mTOR qqnorm montecarlo.png"))

[skewness(montecarloflat[:, i]) for i ∈ 1:4]
[kurtosis(montecarloflat[:, i]) for i ∈ 1:4]


# anova

montecarloresult = anova(montecarlomeasurements[:, 1, :], factornames = ["Rapamycin"])
differencecontrasts(montecarloresult)
