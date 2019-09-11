# Recreates the analysis from the original data files

using SMLMAssociationAnalysis_NCB

datapath =  raw"C:\Users\nicho\Dropbox (Partners HealthCare)\Data Analysis"
projectdirname = "MEG3 Project"
experimentdirnames = ["2 - U2OS p53 MDM2 STORM",
                      "3 - U2OS p53 MEG3 STORM"]

datadirname = "Data"
outputdirname = "Output"

samplenames = ["A", "B", "C", "D"]

nreplicates = 3
nsamples = 4
ncells = 10

experimentresults = Vector{Vector{Vector{Result}}}[]
for experimentdirname ∈ experimentdirnames
    println("Starting experiment $experimentdirname.")
    experimentpath = joinpath(datapath, projectdirname, experimentdirname, datadirname)
    experimentoutputpath = joinpath(datapath, projectdirname, experimentdirname, outputdirname)
    replicateresults = Vector{Vector{Result}}[]
    for i ∈ 1:nreplicates
        sampleresults = Vector{Result}[]
        println("    Starting replicate $i.")
        replicatepath = joinpath(experimentpath, "Replicate $i")
        for samplename ∈ samplenames
            results = Result[]
            println("        Starting sample $samplename.")
            for j ∈ 1:ncells
                println("            Starting cell $j.")
                cellpath = joinpath(replicatepath, "$samplename $(Printf.@sprintf("%03i", j)).bin.txt")
                localizations = LocalizationMicroscopy.load(cellpath, LocalizationMicroscopy.nikonelementstext)
                # account for variances in data collection
                if experimentdirname == experimentdirnames[2] && samplename ∈ samplenames[3:4]
                    ch1_name = "561"
                else
                    ch1_name = "647"
                end
                ch2_name = "488"
                if experimentdirname == experimentdirnames[2] &&
                    ((i == 3 && samplename ∈ samplenames[1:3]) || (i == 2 && samplename == samplenames[2]))
                    ch1_startframe = 1
                    ch2_startframe = 15001
                else
                    ch1_startframe = 1
                    ch2_startframe = 11001 # try weeding out molecules with only 1 loc? ##############
                end
                ch1_molecules, ch1_localizations = getmolecules(localizations, ch1_name, ch1_startframe, 11000, 100, 10, 34.2, 500, 200)
                ch2_molecules, ch2_localizations = getmolecules(localizations, ch2_name, ch2_startframe, 11000, 100, 10, 34.2, 500, 200)
                ch1_neighbors, ch2_neighbors, distances = exclusivenearestneighbors(ch1_molecules, ch2_molecules)

                percentileranks = montecarloaffinity(ch1_molecules, ch2_molecules, ch1_neighbors, ch2_neighbors, distances, 200, 4)

                if length(distances) == 0
                    mediandistance = NaN
                else
                    mediandistance = median(distances)
                end
                println("                Done: $(length(distances)) neighbors from $(length(ch1_molecules)) and $(length(ch2_molecules)) molecules, $(length(ch1_localizations)) and $(length(ch2_localizations)) localizations; median distance $mediandistance")
                ch1_data = ChannelData(ch1_name, ch1_molecules, ch1_neighbors)
                ch2_data = ChannelData(ch2_name, ch2_molecules, ch2_neighbors)
                result = Result(projectdirname, experimentdirname, i, samplename, j,
                                [ch1_data, ch2_data], distances, mediandistance, percentileranks)
                push!(results, result)
            end
            push!(sampleresults, results)
        end
        push!(replicateresults, sampleresults)
    end
    push!(experimentresults, replicateresults)
end

outputpath = joinpath(experimentoutputpath, "experimentresults1.jld2")
save(outputpath, "experimentresults", experimentresults)

#### Analysis
show("Automatically executing this section isn't ideal, as results may not all print, although graphics will be saved. Intended to go through one-by-one as part of the analysis.")

outputdir = joinpath(rootpath, "original", "output")

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

using SimpleANOVA
levene(mediansflat)

# qqnorm, skewness, kurtosis
using StatsPlots
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

using SimpleANOVA
levene(montecarloflat)

# qqnorm, skewness, kurtosis
using StatsPlots
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

using SimpleANOVA
levene(mediansflat)

# qqnorm, skewness, kurtosis
using StatsPlots
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

using SimpleANOVA
levene(montecarloflat)

# qqnorm, skewness, kurtosis
using StatsPlots
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
