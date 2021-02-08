### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ b11fdebe-5507-11eb-1001-537a20e51d0d
begin
	let
		env = mktempdir()
		import Pkg
		Pkg.activate(env)
		Pkg.Registry.update()
		Pkg.add("FileIO")
		Pkg.add("JLD2")
		Pkg.add("StatsBase")
		Pkg.add("StatsPlots")
		Pkg.add("SimpleANOVA")
		Pkg.add(path="../..")
	end
	using SMLMAssociationAnalysis_NCB
	using FileIO
	using StatsBase
	using Statistics
	using StatsPlots
	using SimpleANOVA
end

# ╔═╡ 56dcb910-5507-11eb-0b0a-a1c2aef89207
md"""
# Analysis of Experiment 2: p53-MEG3
"""

# ╔═╡ 90d2acb0-5507-11eb-39c4-df93c258dd98
md"""Load packages"""

# ╔═╡ 86502ff0-550d-11eb-1c13-553251c061ee
md"""Load results file and compute normalized values"""

# ╔═╡ 5e3db410-550d-11eb-0331-49a447464b07
begin
	### Load saved data

	nreplicates = 3
	nsamples = 4
	ncells = 10

	outputdir = "../../output"
	datapath = joinpath(outputdir, "results.jld2")

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

	()
end

# ╔═╡ 294e6a02-550e-11eb-3450-d9d5989ec08c
md"""
Relevant data are stored in these variables:
- `medianmeasurements` - median distance between the pairs
- `montecarlomeasurements` - fraction of pairs considered associated (p < 0.1, dist < 200 nm)
- `normalizedmontecarlomeasurements` - fraction of pairs considered associated, normalized bassed on simulations of 0% and 100% binding

The data is stored in matrices with these dimensions:
- 1 - cell
- 2 - experimental replicate
- 3 - condition 1 (1 = -Dox, 2 = +Dox)
- 4 - condition 2 (Exp 1: 1 = -Nut, 2 = +Nut; Exp 2: 1 = MEG3, 2 = GAPDH)
- 5 - experiment (1 = p53-Mdm2, 2 = p53-MEG3)
"""

# ╔═╡ 9a94f370-5511-11eb-1f60-7fe83ed11def
md"""
### Planning
Central Limit Theorem - Only 10 measurements per level, but 30 per condition. Might be close enough for a normal distribution.

Group sizes equal

Doxycycline and RNA type are independent.
"""

# ╔═╡ 572ba5f2-5513-11eb-1937-3f03cec8252e
md"""
### Medians
#### Check unusual cases
"""

# ╔═╡ 864b181e-5513-11eb-0304-8983d37a8170
let
	p1 = boxplot(medianmeasurements[:, :, 1, 1, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
	p2 = boxplot(medianmeasurements[:, :, 2, 1, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
	p3 = boxplot(medianmeasurements[:, :, 1, 2, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
	p4 = boxplot(medianmeasurements[:, :, 2, 2, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))

	plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "MDM2-p53")
end

# ╔═╡ 27e24820-5514-11eb-1415-8711a194892d
md"""
Potential outliers in box plots:
"""

# ╔═╡ f0393540-5514-11eb-2f76-3514cc3bc85c
[medianmeasurements[5,2,1,1,2],
 medianmeasurements[[4,7],2,2,1,2],
 medianmeasurements[2,1,1,2,2],
medianmeasurements[6,1,2,2,2]]

# ╔═╡ f02df530-5515-11eb-0d9e-8dbb63ed2023
md"""Evaluating z-scores within each replicate:"""

# ╔═╡ ffa9f6c0-576e-11eb-0f11-3139a34ff859
mediansflat = [medianmeasurements[:, :, 1, 1, 2] medianmeasurements[:, :, 2, 1, 2] medianmeasurements[:, :, 1, 2, 2] medianmeasurements[:, :, 2, 2, 2]]

# ╔═╡ 25253910-5516-11eb-0a7a-1d73ed51d4a1
medianzscores = [zscore(mediansflat[:,i]) for i ∈ axes(mediansflat, 2)]

# ╔═╡ d8c96c5e-5517-11eb-0861-056d5e425f04
md"""
- In [1,1] Replicate 1: All within 1.96
- In [1,1] Replicate 2: All within 1.96
- In [1,1] Replicate 3: All within 1.96
- In [1,2] Replicate 1: 1 value, Cell 5 (10%) just under 2.58. None outside 3.29.
- In [1,2] Replicate 2: 1 value, Cell 4 (10%) outside 1.96. None outside 2.58.
- In [1,2] Replicate 3: All within 1.96
- In [2,1] Replicate 1: 1 value, Cell 2 (10%) outside 1.96 significantly. None outside 2.58.
- In [2,1] Replicate 2: All within 1.96
- In [2,1] Replicate 3: 1 value, Cell 3 (10%) outside 1.96 significantly. None outside 2.58.
- In [2,2] Replicate 1: 1 value, Cell 6 (10%) outside 1.96 significantly. None outside 2.58.
- In [2,2] Replicate 2: 1 value, Cell 7 (10%) outside 1.96. None outside 2.58.
- In [2,2] Replicate 3: 1 value, Cell 8 (10%) outside 1.96. None outside 2.58.
"""

# ╔═╡ 87eb8f70-5518-11eb-1639-f91b4bcb5c9b
md"""
None of these potential outliers are due to data entry error or measurement error. Sampling a different population is possible but unknown. These potential outliers are relatively modest, though two quite severe. Since an ANOVA will be used, must be careful not to delete points from a balanced design. Therefore, no outliers will be removed. The two severe outliers will be Windsorized.
"""

# ╔═╡ 74a25f10-59ab-11eb-3a45-bf3b0554a07a
begin
	medianmeasurementsw = copy(medianmeasurements)
	mediansflatw = copy(mediansflat)
	medianmeasurementsw[5,1,2,1,2] = medianmeasurements[4,1,2,1,2]
	mediansflatw[5,4] = mediansflat[4,4]
	medianmeasurementsw[6,1,2,2,2] = medianmeasurements[8,1,2,2,2]
	mediansflatw[6,10] = mediansflat[8,10]
end

# ╔═╡ a87aac80-5518-11eb-38d4-5fe90140b820
md"""
#### Check assumption of homoscedasticity
"""

# ╔═╡ f6e23cd0-5518-11eb-226a-c5751eb65b34
md"""
Based on the ZResid/ZPred plot, the data appears mildly heteroscedastic; outliers really stand out.

Conducted Levene's test, a 1-way ANOVA of the absolute deviation between each value and the group mean, by replicate.
"""

# ╔═╡ 0ffdd990-5519-11eb-376b-33be0cc966d7
levene(mediansflat)

# ╔═╡ 28a8da30-5519-11eb-2fb0-efe693c6e44c
md"""
By Levene's test, the data is very likely heteroscedastic.
 
Group sizes are equal, but samples are on the smaller side. CIs and significance tests may be affected.
"""

# ╔═╡ 5f8de452-5519-11eb-0a2c-051d3e701dd5
md"""
#### Check assumption of normality
"""

# ╔═╡ c4dbef60-55bd-11eb-1822-49c3473000e8
md"""
There are a couple notable deviations but for the *most* part they appear close to normal. The worst ones had the potential outliers identified previously.
"""

# ╔═╡ cf59a26e-55bd-11eb-1c4b-316a603b9007
medianskewness = [skewness(mediansflat[:,i]) for i ∈ axes(mediansflat, 2)]

# ╔═╡ 135dbe6e-55be-11eb-3a89-dd3d02c3ad4a
mediankurtosis = [kurtosis(mediansflat[:,i]) for i ∈ axes(mediansflat, 2)]

# ╔═╡ 6778a880-55be-11eb-1016-c59335beef0a
md"""
From skewness and kurtosis statistics, all were within 1.96 (5% level).
"""

# ╔═╡ 9bf9aaf0-55be-11eb-17a2-75b3eec550f8
md"""
The primary source of non-normality appear to be outliers, relatively minimal. As a result, I don't intend to do anything special.
"""

# ╔═╡ 265006c0-55c1-11eb-2ba8-3d628d096704
md"""
#### Conduct ANOVA tests
"""

# ╔═╡ c101d900-57ab-11eb-04ab-77cd4e0ff4b3
medianresult = anova(
    permutedims(medianmeasurements[:, :, :, :, 2], (1,3,4,2)),
    [fixed, fixed, subject],
    factornames = ["Doxycycline", "RNA", "Replicate"],
)

# ╔═╡ 19934eb0-57ba-11eb-3f45-11c68057c1f7
md"""
#### Check presence of interactions
"""

# ╔═╡ 25bb5660-57ba-11eb-3999-4bbf47d29d8a
plot(medianresult)

# ╔═╡ 342cfff0-57ba-11eb-07d6-2b90717f7c7b
md"""
There seem to be large interactions, but they aren't significant. Since we expect RNAs to be different anyway, will continue to examine ANOVA within each RNA.
"""

# ╔═╡ fd2d4a40-57ab-11eb-2688-53b009a37680
md"""
#### Interpret main effects
RNA identity is not significant and has an insignificant effect size.

Doxycycline is not significant and has an insignificant effect size. As the RNA should be independent, we'll do another ANOVA within each RNA.
"""

# ╔═╡ 5f1ce3b0-55c1-11eb-3d03-b53516cfc976
medianresultMEG3 = anova(
    permutedims(medianmeasurements[:, :, :, 1, 2], (1, 3, 2)),
    [fixed, subject],
    factornames = ["Doxycycline", "Replicate"],
)

# ╔═╡ 0fddbb50-5772-11eb-01cc-e57dbec14a6b
medianresultGAPDH = anova(
    permutedims(medianmeasurements[:, :, :, 2, 2], (1, 3, 2)),
    [fixed, subject],
    factornames = ["Doxycycline", "Replicate"],
)

# ╔═╡ 2ed8ff90-568e-11eb-10a3-0b97b2425d2f
md"""
#### Interpret main effects of within-RNA ANOVAs

Doxycycline produced an insignificant effect with GAPDH (both in size and probability).

Doxycycline produced a large but not significant effect with MEG3, indicating MEG3 and p53 were brought closer together but not unambiguously. This is probably mostly due to the induced expression of MEG3.
"""

# ╔═╡ 25d0ce1e-5691-11eb-2cc8-49980b6b2161
md"""
### Monte Carlo fraction associated
#### Check unusual cases
"""

# ╔═╡ 8c0ebd10-56a9-11eb-3e8c-95b53bf12061
let
	p1 = boxplot(montecarlomeasurements[:, :, 1, 1, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p2 = boxplot(montecarlomeasurements[:, :, 2, 1, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p3 = boxplot(montecarlomeasurements[:, :, 1, 2, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p4 = boxplot(montecarlomeasurements[:, :, 2, 2, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))

	plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "Mdm2-p53")
end

# ╔═╡ faa812c0-56aa-11eb-3b6b-99fe5f200b67
md"""
Potential outliers in box plots:
"""

# ╔═╡ 012a7c50-56ab-11eb-3d50-c9820925f9ac
[montecarlomeasurements[2,3,1,1,2],
 montecarlomeasurements[[7,9],3,1,2,2]]

# ╔═╡ 1637a040-56ac-11eb-2433-49c37383e005
md"""Evaluating z-scores within each replicate:"""

# ╔═╡ 87ec3520-56ac-11eb-2a1a-2981155d4168
montecarloflat = [montecarlomeasurements[:, :, 1, 1, 2] montecarlomeasurements[:, :, 2, 1, 2] montecarlomeasurements[:, :, 1, 2, 2] montecarlomeasurements[:, :, 2, 2, 2]]

# ╔═╡ 9ce1f5f0-56ac-11eb-22a5-d9b5d6a18978
montecarlozscores = [zscore(montecarloflat[:,i]) for i ∈ axes(montecarloflat, 2)]

# ╔═╡ ba7d2760-56ac-11eb-11ca-d543a628c68b
md"""
- In [1,1] Replicate 1: 1 value, Cell 4 (10%) outside 1.96 barely. None outside 2.58.
- In [1,1] Replicate 2: All within 1.96
- In [1,1] Replicate 3: 1 value, Cell 2 (10%) outside 1.96 significantly. None outside 2.58.
- In [1,2] Replicate 1: 1 value, Cell 4 (10%) outside 1.96. None outside 2.58.
- In [1,2] Replicate 2: All within 1.96
- In [1,2] Replicate 3: All within 1.96
- In [2,1] Replicate 1: All within 1.96
- In [2,1] Replicate 2: All within 1.96
- In [2,1] Replicate 3: All within 1.96
- In [2,2] Replicate 1: All within 1.96
- In [2,2] Replicate 2: All within 1.96
- In [2,2] Replicate 3: 1 value, Cell 8 (10%) outside 1.96. None outside 2.58.
"""

# ╔═╡ c18cc050-56ad-11eb-261c-4177aa31ca7c
md"""
None of these potential outliers are due to data entry error or measurement error. Sampling a different population is possible but unknown. These potential outliers are relatively modest. In addition, since an ANOVA will be used, must be careful not to delete points from a balanced design. Therefore, no outliers will be removed.
"""

# ╔═╡ 2bba9b02-56ae-11eb-2e30-439dc20273f8
md"""
#### Check assumption of homoscedasticity
"""

# ╔═╡ 0f8a00a0-56af-11eb-2dc6-abd2a45298a4
md"""
Based on the ZResid/ZPred plot, the data appears very heteroscedastic.

Conducted Levene's test, a 1-way ANOVA of the absolute deviation between each value and the group mean, by replicate.
"""

# ╔═╡ b1ae0750-56af-11eb-0714-35b78e197b5e
levene(montecarloflat)

# ╔═╡ f9676320-56af-11eb-1481-ef9502edded8
md"""
By Levene's test, the data is very likely heteroscedastic.
 
Group sizes are equal, but samples are on the smaller side. CIs and significance tests may be affected.
"""

# ╔═╡ 4e36f730-56b0-11eb-21e8-0d2e7d0a2562
md"""
#### Check assumption of normality
"""

# ╔═╡ 4423b3d0-56b2-11eb-22be-fd86deabef7f
md"""
There are a couple notable deviations but for the *most* part they appear close to normal.
"""

# ╔═╡ db6922c0-56b2-11eb-2601-5525d3f6bf70
montecarloskewness = [skewness(montecarloflat[:,i]) for i ∈ axes(montecarloflat, 2)]

# ╔═╡ 047df370-56b3-11eb-3036-2393031bd967
montecarlokurtosis = [kurtosis(montecarloflat[:,i]) for i ∈ axes(montecarloflat, 2)]

# ╔═╡ 3316f790-56b3-11eb-3d9c-e9ebf797c923
md"""
From skewness and kurtosis statistics, none was greater in value than 1.96 (5% level).
"""

# ╔═╡ 6f1da2c0-56b3-11eb-2506-5f9a944735b6
md"""
#### Conduct ANOVA tests
"""

# ╔═╡ c29f8040-57a7-11eb-2cb4-7974bcf74e27
montecarloresult = anova(
    permutedims(montecarlomeasurements[:, :, :, :, 2], (1, 3, 4, 2)),
    [fixed, fixed, random],
    factornames = ["Doxycycline", "RNA", "Replicate"],
)

# ╔═╡ 2553a3a0-57a9-11eb-247d-a161d85d4410
md"""
#### Interpret main effects
RNA identity is significant and has a very large effect size.

Doxycycline is not significant and has a small effect size. As the RNA should be independent, we'll do another ANOVA within each RNA.
"""

# ╔═╡ 755d8510-56b3-11eb-3a0c-53b643ba2053
montecarloresultMEG3 = anova(
    permutedims(montecarlomeasurements[:, :, :, 1, 2], (1, 3, 2)),
    [fixed, subject],
    factornames = ["Doxycycline", "Replicate"],
)

# ╔═╡ 28212d70-57a7-11eb-1c94-3fa38021ba8e
montecarloresultGAPDH = anova(
    permutedims(montecarlomeasurements[:, :, :, 2, 2], (1, 3, 2)),
    [fixed, subject],
    factornames = ["Doxycycline", "Replicate"],
)

# ╔═╡ 248d1560-5753-11eb-2156-89cfb6742205
md"""
#### Interpret main effects of within-RNA ANOVAs
Doxycycline had a medium effect on MEG3 and p53, although not significant. Doxycyline increased fraction bound (7.18% to 10.3%).

Doxycyline had no effect on GAPDH and p53 and was highly nonsignificant (4.30% to 4.22%).
"""

# ╔═╡ 78883b10-5756-11eb-32cf-b95fdeebda79
md"""
### Normalized Monte Carlo fraction associated
This approach was abandoned because we determined it was over-normalizing the data, but still included here for completeness.

#### Check unusual cases
"""

# ╔═╡ df7dce20-5756-11eb-3796-037e8e3c49f9
let
	p1 = boxplot(normalizedmontecarlomeasurements[:, :, 1, 1, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p2 = boxplot(normalizedmontecarlomeasurements[:, :, 2, 1, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p3 = boxplot(normalizedmontecarlomeasurements[:, :, 1, 2, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p4 = boxplot(normalizedmontecarlomeasurements[:, :, 2, 2, 2], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))

	plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "Mdm2-p53")
end

# ╔═╡ faf5aa60-5756-11eb-3b5b-f1ad5e004af2
md"""
Potential outliers in box plots:
"""

# ╔═╡ ff845a40-5756-11eb-2ce8-d5ba8458e2ae
[normalizedmontecarlomeasurements[4,1,1,1,2],
 normalizedmontecarlomeasurements[9,2,1,1,2],
 normalizedmontecarlomeasurements[5,1,1,2,2],
 normalizedmontecarlomeasurements[7,3,1,2,2]]

# ╔═╡ d2282670-5757-11eb-0bd9-b1aaf0b37d70
md"""Evaluating z-scores within each replicate:"""

# ╔═╡ d99fec80-5757-11eb-39cd-4b0febbf61e3
normmontecarloflat = [normalizedmontecarlomeasurements[:, :, 1, 1, 2] normalizedmontecarlomeasurements[:, :, 2, 1, 2] normalizedmontecarlomeasurements[:, :, 1, 2, 2] normalizedmontecarlomeasurements[:, :, 2, 2, 2]]

# ╔═╡ edee22b2-5757-11eb-317a-d3e4d31486cd
normmontecarlozscores = [zscore(normmontecarloflat[:,i]) for i ∈ axes(normmontecarloflat, 2)]

# ╔═╡ a05efb90-5758-11eb-10ca-fb9fb8e05821
md"""
- In [1,1] Replicate 1: 1 value, Cell 4 (10%) outside 1.96. None outside 2.58.
- In [1,1] Replicate 2: 1 value, Cell 9 (10%) outside 1.96. None outside 2.58.
- In [1,1] Replicate 3: 1 value, Cell 2 (10%) outside 1.96. None outside 2.58.
- In [1,2] Replicate 1: 1 value, Cell 5 (10%) outside 1.96. None outside 2.58.
- In [1,2] Replicate 2: All within 1.96
- In [1,2] Replicate 3: All within 1.96
- In [2,1] Replicate 1: 1 value, Cell 6 (10%) outside 1.96 barely. None outside 2.58.
- In [2,1] Replicate 2: All within 1.96
- In [2,1] Replicate 3: 1 value, Cell 7 (10%) outside 1.96. None outside 2.58.
- In [2,2] Replicate 1: All within 1.96
- In [2,2] Replicate 2: All within 1.96
- In [2,2] Replicate 3: 1 value, Cell 8 (10%) outside 1.96. None outside 2.58.
"""

# ╔═╡ 7323e450-5759-11eb-1c62-7ba711bf418b
md"""
None of these potential outliers are due to data entry error or measurement error. Sampling a different population is possible but unknown. These potential outliers are relatively modest. In addition, since an ANOVA will be used, must be careful not to delete points from a balanced design. Therefore, no outliers will be removed.
"""

# ╔═╡ df4eee90-5759-11eb-02ef-53a0c6f31f07
md"""
#### Check assumption of homoscedasticity
"""

# ╔═╡ 8c2a8790-575b-11eb-34be-99f584826d43
md"""
Based on the ZResid/ZPred plot, the data appears heteroscedastic, although not much relationship between residual and predicted.

Conducted Levene's test, a 1-way ANOVA of the absolute deviation between each value and the group mean, by replicate.
"""

# ╔═╡ 9aade460-575b-11eb-29da-4ff6881ceb60
levene(normmontecarloflat)

# ╔═╡ b8f00430-575b-11eb-0ada-7f0001326feb
md"""
By Levene's test, the data is likely heteroscedastic.
 
Group sizes are equal, but samples are on the smaller side. CIs and significance tests may be affected.
"""

# ╔═╡ 1326cdce-575c-11eb-0f72-cf6eebe27bdf
md"""
#### Check assumption of normality
"""

# ╔═╡ 3c6921b0-575d-11eb-160f-6f06db5965f3
md"""
There are a couple notable deviations but for the *most* part they appear close to normal.
"""

# ╔═╡ 59314bb0-575d-11eb-2b54-e1bf609e6d36
normmontecarloskewness = [skewness(normmontecarloflat[:,i]) for i ∈ axes(normmontecarloflat, 2)]

# ╔═╡ 682c7090-575d-11eb-3b12-d50d450d3fce
normmontecarlokurtosis = [kurtosis(normmontecarloflat[:,i]) for i ∈ axes(normmontecarloflat, 2)]

# ╔═╡ 80f22522-575d-11eb-211b-59fe55838280
md"""
From skewness and kurtosis statistics, none was greater in value than 1.96 (5% level).
"""

# ╔═╡ 91c676d0-575d-11eb-2fd5-97ca14bfdec0
md"""
#### Conduct ANOVA tests
Regular ANOVA
"""

# ╔═╡ b7624d5e-575d-11eb-1221-974537c3b48e
normmontecarloresult = anova(
    permutedims(normalizedmontecarlomeasurements[:, :, :, :, 2], (1, 3, 4, 2)),
    [fixed, fixed, subject],
    factornames = ["Doxycycline", "RNA", "Replicate"],
)

# ╔═╡ 8f0dad80-575f-11eb-2312-adae92d50e83
md"""
### Check presence of interactions
"""

# ╔═╡ a5489880-575f-11eb-1afa-c3bfaf5dccc8
plot(normmontecarloresult)

# ╔═╡ 4b211122-575f-11eb-0aa9-59fc4e9a2df7
md"""
The lines cross, but the F-statistic is not significant. Interpreting simple effects valid. Continuing to look within each RNA.
"""

# ╔═╡ f325b3c0-57ba-11eb-1cae-17d58dd4ec29
normmontecarloresultMEG3 = anova(
    permutedims(normalizedmontecarlomeasurements[:, :, :, 1, 2], (1, 3, 2)),
    [fixed, subject],
    factornames = ["Doxycycline", "Replicate"],
)

# ╔═╡ 0f2bc640-57bb-11eb-2632-e1058fbc9108
normmontecarloresultGAPDH = anova(
    permutedims(normalizedmontecarlomeasurements[:, :, :, 2, 2], (1, 3, 2)),
    [fixed, subject],
    factornames = ["Doxycycline", "Replicate"],
)

# ╔═╡ 408e5d30-5764-11eb-124f-4f2630d1d908
md"""
#### Interpret main effects
Doxycycline had a small but nonsignificant effect with MEG3 (1.55% to 3.72%). Increasing MEG3 concentration may slightly increase the binding of MEG3 and p53, but evidence is weak.

Doxycyline had no significant effect with GAPDH (0.678% to 1.19%). Doxycycline shouldn't do anything to GAPDH, so this is expected.
"""

# ╔═╡ 9c92b290-56ae-11eb-2595-85845fe95f0e
md"""
## Functions
"""

# ╔═╡ 2623d480-56af-11eb-2ab8-0d0afd8c45dc
function zresid_zpred(data)
	z = zscore(data)
	zpred = repeat(mean(z, dims = 1), 10)
	zresid = z .- zpred
	zresid, zpred
end

# ╔═╡ a655a120-56ae-11eb-30be-459d7c3067dc
function zresid_zpred_plot(data)
	zresid, zpred = zresid_zpred(data)
	scatter(
		zresid,
		zpred,
		xaxis = ("Standardized Residual (ZResid)"),
		yaxis = ("Standardized Predicted Value (ZPred)"),
		legend = :none,
	)
end

# ╔═╡ c8f83d60-5518-11eb-2643-a30df2dd9fb9
zresid_zpred_plot(mediansflat)

# ╔═╡ 3a3fa580-56ae-11eb-1f19-f3c06ef2d9d7
zresid_zpred_plot(montecarloflat)

# ╔═╡ 59522180-575a-11eb-39c4-f9eaa27f7901
zresid_zpred_plot(normmontecarloflat)

# ╔═╡ 7922d1d0-56b0-11eb-261c-5f87ca587394
function qqnormplot(data)
	zresid, _ = zresid_zpred(data)
	p = [qqnorm(zresid[:, i]) for i ∈ 1:12]
	plot(p..., layout = grid(4, 3), legend = :none, plot_title = "MDM2-p53 qqnorm")
end

# ╔═╡ 6f727ed0-5519-11eb-1462-01c83aa975b6
qqnormplot(mediansflat)

# ╔═╡ 5ef3f050-56b0-11eb-31f9-419a4161b67e
qqnormplot(montecarloflat)

# ╔═╡ 1823a8d0-575c-11eb-0fa1-893a914f69f7
qqnormplot(normmontecarloflat)

# ╔═╡ Cell order:
# ╟─56dcb910-5507-11eb-0b0a-a1c2aef89207
# ╟─90d2acb0-5507-11eb-39c4-df93c258dd98
# ╠═b11fdebe-5507-11eb-1001-537a20e51d0d
# ╟─86502ff0-550d-11eb-1c13-553251c061ee
# ╠═5e3db410-550d-11eb-0331-49a447464b07
# ╟─294e6a02-550e-11eb-3450-d9d5989ec08c
# ╟─9a94f370-5511-11eb-1f60-7fe83ed11def
# ╟─572ba5f2-5513-11eb-1937-3f03cec8252e
# ╟─864b181e-5513-11eb-0304-8983d37a8170
# ╟─27e24820-5514-11eb-1415-8711a194892d
# ╟─f0393540-5514-11eb-2f76-3514cc3bc85c
# ╟─f02df530-5515-11eb-0d9e-8dbb63ed2023
# ╟─ffa9f6c0-576e-11eb-0f11-3139a34ff859
# ╟─25253910-5516-11eb-0a7a-1d73ed51d4a1
# ╟─d8c96c5e-5517-11eb-0861-056d5e425f04
# ╟─87eb8f70-5518-11eb-1639-f91b4bcb5c9b
# ╠═74a25f10-59ab-11eb-3a45-bf3b0554a07a
# ╟─a87aac80-5518-11eb-38d4-5fe90140b820
# ╠═c8f83d60-5518-11eb-2643-a30df2dd9fb9
# ╟─f6e23cd0-5518-11eb-226a-c5751eb65b34
# ╟─0ffdd990-5519-11eb-376b-33be0cc966d7
# ╟─28a8da30-5519-11eb-2fb0-efe693c6e44c
# ╟─5f8de452-5519-11eb-0a2c-051d3e701dd5
# ╟─6f727ed0-5519-11eb-1462-01c83aa975b6
# ╟─c4dbef60-55bd-11eb-1822-49c3473000e8
# ╟─cf59a26e-55bd-11eb-1c4b-316a603b9007
# ╟─135dbe6e-55be-11eb-3a89-dd3d02c3ad4a
# ╟─6778a880-55be-11eb-1016-c59335beef0a
# ╟─9bf9aaf0-55be-11eb-17a2-75b3eec550f8
# ╟─265006c0-55c1-11eb-2ba8-3d628d096704
# ╠═c101d900-57ab-11eb-04ab-77cd4e0ff4b3
# ╟─19934eb0-57ba-11eb-3f45-11c68057c1f7
# ╟─25bb5660-57ba-11eb-3999-4bbf47d29d8a
# ╟─342cfff0-57ba-11eb-07d6-2b90717f7c7b
# ╟─fd2d4a40-57ab-11eb-2688-53b009a37680
# ╠═5f1ce3b0-55c1-11eb-3d03-b53516cfc976
# ╠═0fddbb50-5772-11eb-01cc-e57dbec14a6b
# ╟─2ed8ff90-568e-11eb-10a3-0b97b2425d2f
# ╟─25d0ce1e-5691-11eb-2cc8-49980b6b2161
# ╟─8c0ebd10-56a9-11eb-3e8c-95b53bf12061
# ╟─faa812c0-56aa-11eb-3b6b-99fe5f200b67
# ╟─012a7c50-56ab-11eb-3d50-c9820925f9ac
# ╟─1637a040-56ac-11eb-2433-49c37383e005
# ╟─87ec3520-56ac-11eb-2a1a-2981155d4168
# ╟─9ce1f5f0-56ac-11eb-22a5-d9b5d6a18978
# ╟─ba7d2760-56ac-11eb-11ca-d543a628c68b
# ╟─c18cc050-56ad-11eb-261c-4177aa31ca7c
# ╟─2bba9b02-56ae-11eb-2e30-439dc20273f8
# ╟─3a3fa580-56ae-11eb-1f19-f3c06ef2d9d7
# ╟─0f8a00a0-56af-11eb-2dc6-abd2a45298a4
# ╟─b1ae0750-56af-11eb-0714-35b78e197b5e
# ╟─f9676320-56af-11eb-1481-ef9502edded8
# ╟─4e36f730-56b0-11eb-21e8-0d2e7d0a2562
# ╟─5ef3f050-56b0-11eb-31f9-419a4161b67e
# ╟─4423b3d0-56b2-11eb-22be-fd86deabef7f
# ╟─db6922c0-56b2-11eb-2601-5525d3f6bf70
# ╟─047df370-56b3-11eb-3036-2393031bd967
# ╟─3316f790-56b3-11eb-3d9c-e9ebf797c923
# ╟─6f1da2c0-56b3-11eb-2506-5f9a944735b6
# ╠═c29f8040-57a7-11eb-2cb4-7974bcf74e27
# ╟─2553a3a0-57a9-11eb-247d-a161d85d4410
# ╠═755d8510-56b3-11eb-3a0c-53b643ba2053
# ╠═28212d70-57a7-11eb-1c94-3fa38021ba8e
# ╟─248d1560-5753-11eb-2156-89cfb6742205
# ╟─78883b10-5756-11eb-32cf-b95fdeebda79
# ╟─df7dce20-5756-11eb-3796-037e8e3c49f9
# ╟─faf5aa60-5756-11eb-3b5b-f1ad5e004af2
# ╠═ff845a40-5756-11eb-2ce8-d5ba8458e2ae
# ╟─d2282670-5757-11eb-0bd9-b1aaf0b37d70
# ╟─d99fec80-5757-11eb-39cd-4b0febbf61e3
# ╟─edee22b2-5757-11eb-317a-d3e4d31486cd
# ╟─a05efb90-5758-11eb-10ca-fb9fb8e05821
# ╟─7323e450-5759-11eb-1c62-7ba711bf418b
# ╟─df4eee90-5759-11eb-02ef-53a0c6f31f07
# ╟─59522180-575a-11eb-39c4-f9eaa27f7901
# ╟─8c2a8790-575b-11eb-34be-99f584826d43
# ╟─9aade460-575b-11eb-29da-4ff6881ceb60
# ╟─b8f00430-575b-11eb-0ada-7f0001326feb
# ╟─1326cdce-575c-11eb-0f72-cf6eebe27bdf
# ╟─1823a8d0-575c-11eb-0fa1-893a914f69f7
# ╟─3c6921b0-575d-11eb-160f-6f06db5965f3
# ╟─59314bb0-575d-11eb-2b54-e1bf609e6d36
# ╟─682c7090-575d-11eb-3b12-d50d450d3fce
# ╟─80f22522-575d-11eb-211b-59fe55838280
# ╟─91c676d0-575d-11eb-2fd5-97ca14bfdec0
# ╠═b7624d5e-575d-11eb-1221-974537c3b48e
# ╟─8f0dad80-575f-11eb-2312-adae92d50e83
# ╟─a5489880-575f-11eb-1afa-c3bfaf5dccc8
# ╟─4b211122-575f-11eb-0aa9-59fc4e9a2df7
# ╠═f325b3c0-57ba-11eb-1cae-17d58dd4ec29
# ╠═0f2bc640-57bb-11eb-2632-e1058fbc9108
# ╟─408e5d30-5764-11eb-124f-4f2630d1d908
# ╟─9c92b290-56ae-11eb-2595-85845fe95f0e
# ╟─2623d480-56af-11eb-2ab8-0d0afd8c45dc
# ╟─a655a120-56ae-11eb-30be-459d7c3067dc
# ╟─7922d1d0-56b0-11eb-261c-5f87ca587394
