# Recreates the analysis from the original data files

using SMLMAssociationAnalysis_NCB

datapath =  raw"C:\Users\nicho\Dropbox (Partners HealthCare)\Data Analysis"
projectdirname = "MEG3 Project"
experimentdirnames = ["7 - U2OS FKBP12 mTOR STORM"]

datadirname = "Data"
outputdirname = "Output"

samplenames = ["E", "F", "G", "H"]

nreplicates = 1
nsamples = 4
ncells = 10

experimentresults = Vector{Vector{Vector{StormData.Result}}}[]
for experimentdirname ∈ experimentdirnames
    println("Starting experiment $experimentdirname.")
    experimentpath = joinpath(datapath, projectdirname, experimentdirname, datadirname)
    experimentoutputpath = joinpath(datapath, projectdirname, experimentdirname, outputdirname)
    replicateresults = Vector{Vector{StormData.Result}}[]
    for i ∈ 1:nreplicates
        sampleresults = Vector{StormData.Result}[]
        println("    Starting replicate $i.")
        replicatepath = joinpath(experimentpath, "Replicate $i")
        for samplename ∈ samplenames
            results = StormData.Result[]
            println("        Starting sample $samplename.")
            for j ∈ 1:ncells
                println("            Starting cell $j.")
                if samplename ∈ samplenames[1:2]
                    # E and F numbers are between 11 and 20 instead of 1 and 10.
                    j = j + 10
                end

                cellpath = joinpath(replicatepath, "$samplename $(Printf.@sprintf("%03i", j)).bin.txt")
                localizations = LocalizationMicroscopy.load(cellpath, LocalizationMicroscopy.nikonelementstext)

                ch1_name = "647"
                ch2_name = "488"
                ch1_startframe = 1
                ch2_startframe = 11001

                ch1_molecules, ch1_localizations = StormAnalysis.getmolecules(localizations, ch1_name, ch1_startframe, 11000, 100, 10, 34.2, 500, 200)
                ch2_molecules, ch2_localizations = StormAnalysis.getmolecules(localizations, ch2_name, ch2_startframe, 11000, 100, 10, 34.2, 500, 200)
                ch1_neighbors, ch2_neighbors, distances = StormAssociation.exclusivenearestneighbors(ch1_molecules, ch2_molecules)

                percentileranks = StormAssociation.montecarloaffinity(ch1_molecules, ch2_molecules, ch1_neighbors, ch2_neighbors, distances, 200, 4)

                if length(distances) == 0
                    mediandistance = NaN
                else
                    mediandistance = median(distances)
                end
                println("                Done: $(length(distances)) neighbors from $(length(ch1_molecules)) and $(length(ch2_molecules)) molecules, $(length(ch1_localizations)) and $(length(ch2_localizations)) localizations; median distance $mediandistance")
                ch1_data = StormData.ChannelData(ch1_name, ch1_molecules, ch1_neighbors)
                ch2_data = StormData.ChannelData(ch2_name, ch2_molecules, ch2_neighbors)
                result = StormData.Result(projectdirname, experimentdirname, i, samplename, j,
                                [ch1_data, ch2_data], distances, mediandistance, percentileranks)
                push!(results, result)
            end
            push!(sampleresults, results)
        end
        push!(replicateresults, sampleresults)
    end
    push!(experimentresults, replicateresults)

    outputpath = joinpath(experimentoutputpath, "results.jld2")
    save(outputpath, "replicateresults", experimentresults)
end



# draw distance-probability plots
for k ∈ 1
    experimentdirname = experimentdirnames[k]
    experimentpath = joinpath(datapath, projectdirname, experimentdirname, datadirname)
    experimentoutputpath = joinpath(datapath, projectdirname, experimentdirname, outputdirname)
    for i ∈ 1:nsamples
        for j ∈ 1:nreplicates
            for m ∈ 1:ncells
                result = experimentresults[k][j][i][m]
                p = StormPlots.distanceprobabilityplot(result)
                Plots.savefig(joinpath(experimentoutputpath, "distprob_$(k)_$(j)_$(i)_$(m).png"))
            end
        end
    end
end

medianmeasurements = Array{Float64, 3}(undef, ncells, nreplicates, 4)
montecarlomeasurements = Array{Float64, 3}(undef, ncells, nreplicates, 4)

for k ∈ 1
    for i ∈ 1:nsamples
        samplemedianresults = Array{Float64, 2}(undef, ncells, nreplicates)
        samplelessthan10results = Array{Float64, 2}(undef, ncells, nreplicates)
        for j ∈ 1:nreplicates
            replicateresults = experimentresults[k][j][i]
            percentilerankslessthan10replicate = [x.percentileranks .< 0.1 for x ∈ replicateresults]
            lessthan10 = count.(percentilerankslessthan10replicate) ./ length.(percentilerankslessthan10replicate)
            mediandistances = map(x -> x.mediandistance, replicateresults)
            samplemedianresults[:,j] = mediandistances
            samplelessthan10results[:,j] = lessthan10
        end
        medianmeasurements[:,:,i] = samplemedianresults
        montecarlomeasurements[:,:,i] = samplelessthan10results
    end
end

medianmeasurements_temp = medianmeasurements[:,:,2]
medianmeasurements[:,:,2] = medianmeasurements[:,:,3]
medianmeasurements[:,:,3] = medianmeasurements_temp

montecarlomeasurements_temp = montecarlomeasurements[:,:,2]
montecarlomeasurements[:,:,2] = montecarlomeasurements[:,:,3]
montecarlomeasurements[:,:,3] = montecarlomeasurements_temp

#### Analysis
show("Automatically executing this section isn't ideal, as results may not all print, although graphics will be saved. Intended to go through one-by-one as part of the analysis.")

outputdir = joinpath(rootpath, "original", "output")

# Check for unusual cases
p1 = boxplot(medianmeasurements[:,:,1], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))
p2 = boxplot(medianmeasurements[:,:,2], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))
p3 = boxplot(medianmeasurements[:,:,3], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))
p4 = boxplot(medianmeasurements[:,:,4], xaxis=("Replicates", [1, 2, 3]), yaxis=("Median distance (nm)"))

plot(p1, p2, p3, p4, layout=grid(2,2), legend=:none, plot_title="FKBP12-mTOR")
savefig(joinpath(outputdir, "FKBP12-mTOR boxplots-afterwindsorizing.png"))

# Windsorize an extreme value
medianmeasurements[4,1,1] = medianmeasurements[5,1,1]


#ZResid/ZPred plot and Levene's test
mediansflat = [medianmeasurements[:,:,1] medianmeasurements[:,:,2] medianmeasurements[:,:,3] medianmeasurements[:,:,4]]
z = zscore(mediansflat)
zpred = repeat(mean(z, dims = 1), 10)
zresid = z .- zpred
scatter(zresid,zpred, xaxis=("Standardized Residual (ZResid)"), yaxis=("Standardized Predicted Value (ZPred)"), legend=:none)
savefig(joinpath(outputdir, "FKBP12-mTOR zresid-zpred.png"))

using SimpleANOVA
levene(mediansflat)

# qqnorm, skewness, kurtosis
using StatsPlots
p = [qqnorm(zresid[:,i]) for i ∈ 1:4]
plot(p..., layout=grid(4,1), legend=:none, plot_title="MEG3-p53 qqnorm")
savefig(joinpath(outputdir, "FKBP12-mTOR qqnorm.png"))

[skewness(mediansflat[:,i]) for i ∈ 1:4]
[kurtosis(mediansflat[:,i]) for i ∈ 1:4]


# anova

medianresult = anova(medianmeasurements[:,1,:], factornames = ["Rapamycin"])


# contrasts
contrastcoeffs_0vs10_20_40 = [-3; 1; 1; 1]
contrast0vs10_20_40 = sum(contrastcoeffs_0vs10_20_40 .* medianresult.cellmeans) / sqrt(medianresult.effects[end].ms * sum(contrastcoeffs_0vs10_20_40 .^ 2 ./ 10))
r0vs10_20_40 = sqrt(contrast0vs10_20_40 ^ 2 / (contrast0vs10_20_40 ^ 2 + medianresult.effects[end].df))

contrastcoeffs_10vs20_40 = [-2; 1; 1]
contrast10vs20_40 = sum(contrastcoeffs_10vs20_40 .* medianresult.cellmeans[2:4]) / sqrt(medianresult.effects[end].ms * sum(contrastcoeffs_10vs20_40 .^ 2 ./ 10))
r10vs20_40 = sqrt(contrast10vs20_40 ^ 2 / (contrast10vs20_40 ^ 2 + medianresult.effects[end].df))

contrastcoeffs_20vs40 = [-1; 1]
contrast20vs40 = sum(contrastcoeffs_20vs40 .* medianresult.cellmeans[3:4]) / sqrt(medianresult.effects[end].ms * sum(contrastcoeffs_20vs40 .^ 2 ./ 10))
r20vs40 = sqrt(contrast20vs40 ^ 2 / (contrast20vs40 ^ 2 + medianresult.effects[end].df))

using Distributions
dist = TDist(medianresult.effects[end].df)
pvalue(dist, x) = min(2 * min(cdf(dist, x), ccdf(dist, x)), 1.0)
p0vs10_20_40 = pvalue(dist, contrast0vs10_20_40)
p10vs20_40 = pvalue(dist, contrast10vs20_40)
p20vs40 = pvalue(dist, contrast20vs40)

### Monte Carlo Exp 3

# Check for unusual cases
p1 = boxplot(montecarlomeasurements[:,:,1], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))
p2 = boxplot(montecarlomeasurements[:,:,2], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))
p3 = boxplot(montecarlomeasurements[:,:,3], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))
p4 = boxplot(montecarlomeasurements[:,:,4], xaxis=("Replicates", [1, 2, 3]), yaxis=("Fraction bound"))

plot(p1, p2, p3, p4, layout=grid(2,2), legend=:none, plot_title="FKBP12-mTOR")
savefig(joinpath(outputdir, "FKBP12-mTOR boxplots montecarlo.png"))

#ZResid/ZPred plot and Levene's test
montecarloflat = [montecarlomeasurements[:,:,1] montecarlomeasurements[:,:,2] montecarlomeasurements[:,:,3] montecarlomeasurements[:,:,4]]
z = zscore(montecarloflat)
zpred = repeat(mean(z, dims = 1), 10)
zresid = z .- zpred
scatter(zresid,zpred, xaxis=("Standardized Residual (ZResid)"), yaxis=("Standardized Predicted Value (ZPred)"), legend=:none)
savefig(joinpath(outputdir, "FKBP12-mTOR zresid-zpred montecarlo.png"))

using SimpleANOVA
levene(montecarloflat)

# qqnorm, skewness, kurtosis
using StatsPlots
p = [qqnorm(zresid[:,i]) for i ∈ 1:4]
plot(p..., layout=grid(4,1), legend=:none, plot_title="FKBP12-mTOR qqnorm montecarlo")
savefig(joinpath(outputdir, "FKBP12-mTOR qqnorm montecarlo.png"))

[skewness(montecarloflat[:,i]) for i ∈ 1:4]
[kurtosis(montecarloflat[:,i]) for i ∈ 1:4]


# anova

montecarloresult = anova(montecarlomeasurements[:,1,:], factornames = ["Rapamycin"])

contrastcoeffs_0vs10_20_40 = [-3; 1; 1; 1]
contrast0vs10_20_40 = sum(contrastcoeffs_0vs10_20_40 .* montecarloresult.cellmeans) / sqrt(montecarloresult.effects[end].ms * sum(contrastcoeffs_0vs10_20_40 .^ 2 ./ 10))
r0vs10_20_40 = sqrt(contrast0vs10_20_40 ^ 2 / (contrast0vs10_20_40 ^ 2 + montecarloresult.effects[end].df))

contrastcoeffs_10vs20_40 = [-2; 1; 1]
contrast10vs20_40 = sum(contrastcoeffs_10vs20_40 .* montecarloresult.cellmeans[2:4]) / sqrt(montecarloresult.effects[end].ms * sum(contrastcoeffs_10vs20_40 .^ 2 ./ 10))
r10vs20_40 = sqrt(contrast10vs20_40 ^ 2 / (contrast10vs20_40 ^ 2 + montecarloresult.effects[end].df))

contrastcoeffs_20vs40 = [-1; 1]
contrast20vs40 = sum(contrastcoeffs_20vs40 .* montecarloresult.cellmeans[3:4]) / sqrt(montecarloresult.effects[end].ms * sum(contrastcoeffs_20vs40 .^ 2 ./ 10))
r20vs40 = sqrt(contrast20vs40 ^ 2 / (contrast20vs40 ^ 2 + montecarloresult.effects[end].df))

using Distributions
dist = TDist(montecarloresult.effects[end].df)
pvalue(dist, x) = min(2 * min(cdf(dist, x), ccdf(dist, x)), 1.0)
p0vs10_20_40 = pvalue(dist, contrast0vs10_20_40)
p10vs20_40 = pvalue(dist, contrast10vs20_40)
p20vs40 = pvalue(dist, contrast20vs40)
