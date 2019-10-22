# Recreates the analysis from the original data files

using SMLMAssociationAnalysis_NCB
using FileIO
using StatsBase
using StatsPlots

### Load saved data

rootpath =  raw"C:\Users\nicho\Dropbox (Partners HealthCare)\Data Analysis"

samplenames = ["A", "B", "C", "D"]

nreplicates = 3
nsamples = 4
ncells = 10

outputdir = joinpath(rootpath, "SMLMAssociationAnalysis_NCB.jl", "original", "output")
datapath = joinpath(outputdir, "results.jld2")

experimentresults = FileIO.load(datapath)["experimentresults"]

medianmeasurements = Array{Float64, 4}(undef, ncells, nreplicates, 4, 2)
montecarlomeasurements = Array{Float64, 4}(undef, ncells, nreplicates, 4, 2)
localization1counts = Array{Float64, 4}(undef, ncells, nreplicates, 4, 2)
molecule1counts = Array{Float64, 4}(undef, ncells, nreplicates, 4, 2)
localization2counts = Array{Float64, 4}(undef, ncells, nreplicates, 4, 2)
molecule2counts = Array{Float64, 4}(undef, ncells, nreplicates, 4, 2)

for k ∈ 1:2
    for i ∈ 1:nsamples
        samplemedianresults = Array{Float64, 2}(undef, ncells, nreplicates)
        samplelessthan10results = Array{Float64, 2}(undef, ncells, nreplicates)
        samplelocalization1counts = Array{Float64, 2}(undef, ncells, nreplicates)
        samplemolecule1counts = Array{Float64, 2}(undef, ncells, nreplicates)
        samplelocalization2counts = Array{Float64, 2}(undef, ncells, nreplicates)
        samplemolecule2counts = Array{Float64, 2}(undef, ncells, nreplicates)
        for j ∈ 1:nreplicates
            replicateresults = experimentresults[k][j][i]
            lessthanlimitreplicate = [(x.distances .< 200) .& (x.percentileranks .< 0.1) for x ∈ replicateresults]
            lessthan10 = count.(lessthanlimitreplicate) ./ length.(lessthanlimitreplicate)
            mediandistances = map(x -> x.mediandistance, replicateresults)
            samplemedianresults[:,j] = mediandistances
            samplelessthan10results[:,j] = lessthan10
            samplelocalization1counts[:,j] = [map(x -> x.group.localizations |> length, y.channels[1].molecules) |> sum for y ∈ replicateresults]
            samplelocalization2counts[:,j] = [map(x -> x.group.localizations |> length, y.channels[2].molecules) |> sum for y ∈ replicateresults]
            samplemolecule1counts[:,j] = map(x -> x.channels[1].molecules |> length, replicateresults)
            samplemolecule2counts[:,j] = map(x -> x.channels[2].molecules |> length, replicateresults)
        end
        medianmeasurements[:,:,i,k] = samplemedianresults
        montecarlomeasurements[:,:,i,k] = samplelessthan10results
        localization1counts[:,:,i,k] = samplelocalization1counts
        molecule1counts[:,:,i,k] = samplemolecule1counts
        localization2counts[:,:,i,k] = samplelocalization2counts
        molecule2counts[:,:,i,k] = samplemolecule2counts
    end
end

medianmeasurements = cat(cat(medianmeasurements[:,:,1:2,1], medianmeasurements[:,:,3:4,1], dims = 4),
                         cat(medianmeasurements[:,:,1:2,2], medianmeasurements[:,:,3:4,2], dims = 4), dims = 5)
montecarlomeasurements = cat(cat(montecarlomeasurements[:,:,1:2,1], montecarlomeasurements[:,:,3:4,1], dims = 4),
                             cat(montecarlomeasurements[:,:,1:2,2], montecarlomeasurements[:,:,3:4,2], dims = 4), dims = 5)
localization1counts = cat(cat(localization1counts[:,:,1:2,1], localization1counts[:,:,3:4,1], dims = 4),
                      cat(localization1counts[:,:,1:2,2], localization1counts[:,:,3:4,2], dims = 4), dims = 5)
molecule1counts = cat(cat(molecule1counts[:,:,1:2,1], molecule1counts[:,:,3:4,1], dims = 4),
                          cat(molecule1counts[:,:,1:2,2], molecule1counts[:,:,3:4,2], dims = 4), dims = 5)
localization2counts = cat(cat(localization2counts[:,:,1:2,1], localization2counts[:,:,3:4,1], dims = 4),
                    cat(localization2counts[:,:,1:2,2], localization2counts[:,:,3:4,2], dims = 4), dims = 5)
molecule2counts = cat(cat(molecule2counts[:,:,1:2,1], molecule2counts[:,:,3:4,1], dims = 4),
                        cat(molecule2counts[:,:,1:2,2], molecule2counts[:,:,3:4,2], dims = 4), dims = 5)

#### Analysis
show("Automatically executing this section isn't ideal, as results may not all print, although graphics will be saved. Intended to go through one-by-one as part of the analysis.")

### Median Exp 2

# Check for unusual cases
p1 = boxplot(medianmeasurements[:,:,1,1,1], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))
p2 = boxplot(medianmeasurements[:,:,2,1,1], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))
p3 = boxplot(medianmeasurements[:,:,1,2,1], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))
p4 = boxplot(medianmeasurements[:,:,2,2,1], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))

plot(p1, p2, p3, p4, layout=grid(2,2), legend=:none, plot_title="MDM2-p53")
savefig(joinpath(outputdir, "Mdm2-p53 boxplots.png"))


#ZResid/ZPred plot and Levene's test
mediansflat = [medianmeasurements[:,:,1,1,1] medianmeasurements[:,:,2,1,1] medianmeasurements[:,:,1,2,1] medianmeasurements[:,:,2,2,1]]
z = zscore(mediansflat)
zpred = repeat(mean(z, dims = 1), 10)
zresid = z .- zpred
scatter(zresid,zpred, xaxis=("Standardized Residual (ZResid)"), yaxis=("Standardized Predicted Value (ZPred)"), legend=:none)
savefig(joinpath(outputdir, "Mdm2-p53 zresid-zpred.png"))

levene(mediansflat)

# qqnorm, skewness, kurtosis
p = [qqnorm(zresid[:,i]) for i ∈ 1:12]
plot(p..., layout=grid(4,3), legend=:none, plot_title="MDM2-p53 qqnorm")
savefig(joinpath(outputdir, "Mdm2-p53 qqnorm.png"))

[skewness(mediansflat[:,i]) for i ∈ 1:12]
[kurtosis(mediansflat[:,i]) for i ∈ 1:12]


# anova

medianresult = anova(medianmeasurements[:,:,:,:,1], [nested], factornames = ["Replicate", "Doxycycline", "Nutlin-3a"])

plot(medianresult)
savefig(joinpath(outputdir, "Mdm2-p53 interactionplot.png"))




### Monte Carlo Exp 2

# Check for unusual cases
p1 = boxplot(montecarlomeasurements[:,:,1,1,1], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))
p2 = boxplot(montecarlomeasurements[:,:,2,1,1], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))
p3 = boxplot(montecarlomeasurements[:,:,1,2,1], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))
p4 = boxplot(montecarlomeasurements[:,:,2,2,1], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))

plot(p1, p2, p3, p4, layout=grid(2,2), legend=:none, plot_title="Mdm2-p53")
savefig(joinpath(outputdir, "Mdm2-p53 boxplots montecarlo.png"))


#ZResid/ZPred plot and Levene's test
montecarloflat = [montecarlomeasurements[:,:,1,1,1] montecarlomeasurements[:,:,2,1,1] montecarlomeasurements[:,:,1,2,1] montecarlomeasurements[:,:,2,2,1]]
z = zscore(montecarloflat)
zpred = repeat(mean(z, dims = 1), 10)
zresid = z .- zpred
scatter(zresid,zpred, xaxis=("Standardized Residual (ZResid)"), yaxis=("Standardized Predicted Value (ZPred)"), legend=:none)
savefig(joinpath(outputdir, "Mdm2-p53 zresid-zpred montecarlo.png"))

levene(montecarloflat)

# qqnorm, skewness, kurtosis
p = [qqnorm(zresid[:,i]) for i ∈ 1:12]
plot(p..., layout=grid(4,3), legend=:none, plot_title="MDM2-p53 qqnorm montecarlo")
savefig(joinpath(outputdir, "Mdm2-p53 qqnorm montecarlo.png"))

[skewness(montecarloflat[:,i]) for i ∈ 1:12]
[kurtosis(montecarloflat[:,i]) for i ∈ 1:12]


# anova

montecarloresult = anova(montecarlomeasurements[:,:,:,:,1], [nested], factornames = ["Replicate", "Doxycycline", "Nutlin-3a"])

plot(montecarloresult)
savefig(joinpath(outputdir, "Mdm2-p53 interactionplot montecarlo.png"))


### Median Exp 3

# Check for unusual cases
p1 = boxplot(medianmeasurements[:,:,1,1,2], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))
p2 = boxplot(medianmeasurements[:,:,2,1,2], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))
p3 = boxplot(medianmeasurements[:,:,1,2,2], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))
p4 = boxplot(medianmeasurements[:,:,2,2,2], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))

plot(p1, p2, p3, p4, layout=grid(2,2), legend=:none, plot_title="MEG3-p53")
savefig(joinpath(outputdir, "MEG3-p53 boxplots.png"))


#ZResid/ZPred plot and Levene's test
mediansflat = [medianmeasurements[:,:,1,1,2] medianmeasurements[:,:,2,1,2] medianmeasurements[:,:,1,2,2] medianmeasurements[:,:,2,2,2]]
z = zscore(mediansflat)
zpred = repeat(mean(z, dims = 1), 10)
zresid = z .- zpred
scatter(zresid,zpred, xaxis=("Standardized Residual (ZResid)"), yaxis=("Standardized Predicted Value (ZPred)"), legend=:none)
savefig(joinpath(outputdir, "MEG3-p53 zresid-zpred.png"))

levene(mediansflat)

# qqnorm, skewness, kurtosis
p = [qqnorm(zresid[:,i]) for i ∈ 1:12]
plot(p..., layout=grid(4,3), legend=:none, plot_title="MEG3-p53 qqnorm")
savefig(joinpath(outputdir, "MEG3-p53 qqnorm.png"))

[skewness(mediansflat[:,i]) for i ∈ 1:12]
[kurtosis(mediansflat[:,i]) for i ∈ 1:12]


# anova

medianresult = anova(medianmeasurements[:,:,:,:,2], [nested], factornames = ["Replicate", "Doxycycline", "RNA"])

plot(medianresult)
savefig(joinpath(outputdir, "MEG3-p53 interactionplot.png"))




### Monte Carlo Exp 3

# Check for unusual cases
p1 = boxplot(montecarlomeasurements[:,:,1,1,2], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))
p2 = boxplot(montecarlomeasurements[:,:,2,1,2], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))
p3 = boxplot(montecarlomeasurements[:,:,1,2,2], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))
p4 = boxplot(montecarlomeasurements[:,:,2,2,2], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))

plot(p1, p2, p3, p4, layout=grid(2,2), legend=:none, plot_title="MEG3-p53")
savefig(joinpath(outputdir, "MEG3-p53 boxplots montecarlo.png"))


#ZResid/ZPred plot and Levene's test
montecarloflat = [montecarlomeasurements[:,:,1,1,2] montecarlomeasurements[:,:,2,1,2] montecarlomeasurements[:,:,1,2,2] montecarlomeasurements[:,:,2,2,2]]
z = zscore(montecarloflat)
zpred = repeat(mean(z, dims = 1), 10)
zresid = z .- zpred
scatter(zresid,zpred, xaxis=("Standardized Residual (ZResid)"), yaxis=("Standardized Predicted Value (ZPred)"), legend=:none)
savefig(joinpath(outputdir, "MEG3-p53 zresid-zpred montecarlo.png"))

levene(montecarloflat)

# qqnorm, skewness, kurtosis
p = [qqnorm(zresid[:,i]) for i ∈ 1:12]
plot(p..., layout=grid(4,3), legend=:none, plot_title="MEG3-p53 qqnorm montecarlo")
savefig(joinpath(outputdir, "MEG3-p53 qqnorm montecarlo.png"))

[skewness(montecarloflat[:,i]) for i ∈ 1:12]
[kurtosis(montecarloflat[:,i]) for i ∈ 1:12]


# anova

montecarloresult = anova(montecarlomeasurements[:,:,:,:,2], [nested], factornames = ["Replicate", "Doxycycline", "RNA"])

plot(montecarloresult)
savefig(joinpath(outputdir, "MEG3-p53 interactionplot montecarlo.png"))

montecarloresultMEG3 = anova(montecarlomeasurements[:,:,:,1,2], [nested], factornames = ["Replicate", "Doxycycline"])
montecarloresultGAPDH = anova(montecarlomeasurements[:,:,:,2,2], [nested], factornames = ["Replicate", "Doxycycline"])
