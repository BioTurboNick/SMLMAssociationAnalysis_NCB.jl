# Recreates the analysis from the original data files

using SMLMAssociationAnalysis_NCB
using FileIO
using StatsBase
using StatsPlots

### Load saved data

rootpath = raw"C:\Users\nicho\Dropbox (Partners HealthCare)\Data Analysis"

samplenames = ["E", "F", "G", "H"]

nreplicates = 1
nsamples = 4
ncells = 10

outputdir = joinpath(rootpath, "SMLMAssociationAnalysis_NCB.jl", "original", "control", "output")
datapath = joinpath(outputdir, "results.jld2")

experimentresults = FileIO.load(datapath)["replicateresults"]

medianmeasurements = Array{Float64,3}(undef, ncells, nreplicates, 4)
montecarlomeasurements = Array{Float64,3}(undef, ncells, nreplicates, 4)

for k ∈ 1
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
savefig(joinpath(outputdir, "FKBP12-mTOR boxplots-afterwindsorizing.png"))

# Windsorize an extreme value
medianmeasurements[4, 1, 1] = medianmeasurements[5, 1, 1]


#ZResid/ZPred plot and Levene's test
mediansflat = [medianmeasurements[:, :, 1] medianmeasurements[:, :, 2] medianmeasurements[:, :, 3] medianmeasurements[
    :,
    :,
    4,
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


# contrasts
contrastcoeffs_0vs10_20_40 = [-3; 1; 1; 1]
contrast0vs10_20_40 = sum(contrastcoeffs_0vs10_20_40 .* medianresult.cellmeans) /
                      sqrt(medianresult.effects[end].ms * sum(contrastcoeffs_0vs10_20_40 .^ 2 ./ 10))
r0vs10_20_40 = sqrt(contrast0vs10_20_40^2 / (contrast0vs10_20_40^2 + medianresult.effects[end].df))

contrastcoeffs_10vs20_40 = [-2; 1; 1]
contrast10vs20_40 = sum(contrastcoeffs_10vs20_40 .* medianresult.cellmeans[2:4]) /
                    sqrt(medianresult.effects[end].ms * sum(contrastcoeffs_10vs20_40 .^ 2 ./ 10))
r10vs20_40 = sqrt(contrast10vs20_40^2 / (contrast10vs20_40^2 + medianresult.effects[end].df))

contrastcoeffs_20vs40 = [-1; 1]
contrast20vs40 = sum(contrastcoeffs_20vs40 .* medianresult.cellmeans[3:4]) /
                 sqrt(medianresult.effects[end].ms * sum(contrastcoeffs_20vs40 .^ 2 ./ 10))
r20vs40 = sqrt(contrast20vs40^2 / (contrast20vs40^2 + medianresult.effects[end].df))

using Distributions
dist = TDist(medianresult.effects[end].df)
pvalue(dist, x) = min(2 * min(cdf(dist, x), ccdf(dist, x)), 1.0)
p0vs10_20_40 = pvalue(dist, contrast0vs10_20_40)
p10vs20_40 = pvalue(dist, contrast10vs20_40)
p20vs40 = pvalue(dist, contrast20vs40)

### Monte Carlo Exp 3

# Check for unusual cases
p1 = boxplot(montecarlomeasurements[:, :, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p2 = boxplot(montecarlomeasurements[:, :, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p3 = boxplot(montecarlomeasurements[:, :, 3], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
p4 = boxplot(montecarlomeasurements[:, :, 4], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))

plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "FKBP12-mTOR")
savefig(joinpath(outputdir, "FKBP12-mTOR boxplots montecarlo.png"))

#ZResid/ZPred plot and Levene's test
montecarloflat = [montecarlomeasurements[:, :, 1] montecarlomeasurements[:, :, 2] montecarlomeasurements[:, :, 3] montecarlomeasurements[
    :,
    :,
    4,
]]
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

contrastcoeffs_0vs10_20_40 = [-3; 1; 1; 1]
contrast0vs10_20_40 = sum(contrastcoeffs_0vs10_20_40 .* montecarloresult.cellmeans) /
                      sqrt(montecarloresult.effects[end].ms * sum(contrastcoeffs_0vs10_20_40 .^ 2 ./ 10))
r0vs10_20_40 = sqrt(contrast0vs10_20_40^2 / (contrast0vs10_20_40^2 + montecarloresult.effects[end].df))

contrastcoeffs_10vs20_40 = [-2; 1; 1]
contrast10vs20_40 = sum(contrastcoeffs_10vs20_40 .* montecarloresult.cellmeans[2:4]) /
                    sqrt(montecarloresult.effects[end].ms * sum(contrastcoeffs_10vs20_40 .^ 2 ./ 10))
r10vs20_40 = sqrt(contrast10vs20_40^2 / (contrast10vs20_40^2 + montecarloresult.effects[end].df))

contrastcoeffs_20vs40 = [-1; 1]
contrast20vs40 = sum(contrastcoeffs_20vs40 .* montecarloresult.cellmeans[3:4]) /
                 sqrt(montecarloresult.effects[end].ms * sum(contrastcoeffs_20vs40 .^ 2 ./ 10))
r20vs40 = sqrt(contrast20vs40^2 / (contrast20vs40^2 + montecarloresult.effects[end].df))

using Distributions
dist = TDist(montecarloresult.effects[end].df)
pvalue(dist, x) = min(2 * min(cdf(dist, x), ccdf(dist, x)), 1.0)
p0vs10_20_40 = pvalue(dist, contrast0vs10_20_40)
p10vs20_40 = pvalue(dist, contrast10vs20_40)
p20vs40 = pvalue(dist, contrast20vs40)
