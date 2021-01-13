### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ b11fdebe-5507-11eb-1001-537a20e51d0d
begin
	import Pkg
	Pkg.add(path="../..")
	Pkg.add("FileIO")
	Pkg.add("JLD2")
	Pkg.add("StatsBase")
	Pkg.add("StatsPlots")
	Pkg.add("SimpleANOVA")
	using SMLMAssociationAnalysis_NCB
	using FileIO
	using StatsBase
	using Statistics
	using StatsPlots
	using SimpleANOVA
end

# ╔═╡ 56dcb910-5507-11eb-0b0a-a1c2aef89207
md"""
# Analysis
"""

# ╔═╡ 90d2acb0-5507-11eb-39c4-df93c258dd98
md"""Load packages"""

# ╔═╡ 86502ff0-550d-11eb-1c13-553251c061ee
md"""Load results file and compute normalized values"""

# ╔═╡ 5e3db410-550d-11eb-0331-49a447464b07
begin
	### Load saved data

	samplenames = ["A", "B", "C", "D"]

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
## Experiment 1

### Planning
Central Limit Theorem - Only 10 measurements per level, but 30 per condition. Might be close enough for a normal distribution.

Group sizes equal

Doxycycline and Nutlin treatments should be independent.
"""

# ╔═╡ 572ba5f2-5513-11eb-1937-3f03cec8252e
md"""
### Medians
#### Check unusual cases
"""

# ╔═╡ 864b181e-5513-11eb-0304-8983d37a8170
begin
	p1 = boxplot(medianmeasurements[:, :, 1, 1, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
	p2 = boxplot(medianmeasurements[:, :, 2, 1, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
	p3 = boxplot(medianmeasurements[:, :, 1, 2, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))
	p4 = boxplot(medianmeasurements[:, :, 2, 2, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Median distance (nm)"))

	plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "MDM2-p53")
end

# ╔═╡ 27e24820-5514-11eb-1415-8711a194892d
md"""
Potential outliers in box plots:
"""

# ╔═╡ f0393540-5514-11eb-2f76-3514cc3bc85c
[medianmeasurements[[3,4],2,1,1,1],
 medianmeasurements[8,1,2,1,1],
 medianmeasurements[6,2,2,1,1],
 medianmeasurements[10,3,1,2,1],
 medianmeasurements[10,2,2,2,1],
 medianmeasurements[5,3,2,2,1]]

# ╔═╡ f02df530-5515-11eb-0d9e-8dbb63ed2023
md"""Evaluating z-scores within each replicate:"""

# ╔═╡ 1bf98170-5516-11eb-0c24-77ac43672fba
mediansflat = [medianmeasurements[:, :, 1, 1, 1] medianmeasurements[:, :, 2, 1, 1] medianmeasurements[:, :, 1, 2, 1] medianmeasurements[:, :, 2, 2, 1]]

# ╔═╡ 25253910-5516-11eb-0a7a-1d73ed51d4a1
medianzscores = [zscore(mediansflat[:,i]) for i ∈ axes(mediansflat, 2)]

# ╔═╡ d8c96c5e-5517-11eb-0861-056d5e425f04
md"""
- In [1,1] Replicate 1: All within 1.96
- In [1,1] Replicate 2: All within 1.96
- In [1,1] Replicate 3: All within 1.96
- In [1,2] Replicate 1: 1 value, Cell 8 (10%) outside 1.96, but just barely. None outside 2.58.
- In [1,2] Replicate 2: 1 value, Cell 6 (10%) outside 1.96, by a lot. None outside 2.58.
- In [1,2] Replicate 3: All within 1.96
- In [2,1] Replicate 1: All within 1.96
- In [2,1] Replicate 2: All within 1.96
- In [2,2] Replicate 1: 1 value, Cell 1 (10%) outside 1.96, but just barely. None outside 2.58.
- In [2,2] Replicate 2: 1 value, Cell 10 (10%) outside 2.58 significantly. None outside 3.29.
- In [2,2] Replicate 3: 1 value, Cell 5 (10%) outside 1.96 significantly. None outside 2.58.
"""

# ╔═╡ 87eb8f70-5518-11eb-1639-f91b4bcb5c9b
md"""
None of these potential outliers are due to data entry error or measurement error. Sampling a different population is possible but unknown. These potential outliers are relatively modest. In addition, since an ANOVA will be used, must be careful not to delete points from a balanced design. Therefore, no outliers will be removed. However, I should compare against a robust method to verify results.
"""

# ╔═╡ a87aac80-5518-11eb-38d4-5fe90140b820
md"""
### Check assumption of homoscedasticity
"""

# ╔═╡ c8f83d60-5518-11eb-2643-a30df2dd9fb9
begin
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
end

# ╔═╡ f6e23cd0-5518-11eb-226a-c5751eb65b34
md"""
Based on the ZResid/ZPred plot, the data appears fairly heteroscedastic.

Conducted Levene's test, a 1-way ANOVA of the absolute deviation between each value and the group mean, by replicate.
"""

# ╔═╡ 0ffdd990-5519-11eb-376b-33be0cc966d7
levene(mediansflat)

# ╔═╡ 28a8da30-5519-11eb-2fb0-efe693c6e44c
md"""
By Levene's test, the data is very likely heteroscedastic. Based on median, still significant at p = 0.011.
 
Not doing Hartley's Fmax because it isn't as simple.
 
Group sizes are equal, but samples are on the smaller side. CIs and significance tests may be affected.

I will apply a corrected F test (Welch's), assuming other violations of assumptions do not exist.

Except that performing Welch is impractical for a factorial ANOVA, a robust test is preferred.
"""

# ╔═╡ 5f8de452-5519-11eb-0a2c-051d3e701dd5
md"""
### Check assumption of normality
"""

# ╔═╡ 6f727ed0-5519-11eb-1462-01c83aa975b6
let
	p = [qqnorm(zresid[:, i]) for i ∈ 1:12]
	plot(p..., layout = grid(4, 3), legend = :none, plot_title = "MDM2-p53 qqnorm")
end

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
From skewness and kurtosis statistics, one was greater in absolute value than 1.96 (5% level):
 [2,2] Replicate 2, skewness: 2.10 (moderate); kurtosis: 3.35 (extreme)
"""

# ╔═╡ 9bf9aaf0-55be-11eb-17a2-75b3eec550f8
md"""
Because the sample size is small, used Shapiro-Wilk test (SPSS). [1,1] Replicate 2, [2,1] Replicate 2 are significantly non-normal (p = 0.012, 0.027, respectively); [2,2] Replicate 2 borderline (p = 0.086).

The primary source of non-normality appear to be outliers, relatively minimal. As a result, I don't intend to do anything special other than verify conclusions with a robust test.
"""

# ╔═╡ 265006c0-55c1-11eb-2ba8-3d628d096704
md"""
### Conduct ANOVA tests
Regular ANOVA
"""

# ╔═╡ 5f1ce3b0-55c1-11eb-3d03-b53516cfc976
medianresult = anova(
    medianmeasurements[:, :, :, :, 1],
    [nested],
    factornames = ["Replicate", "Doxycycline", "Nutlin-3a"],
)

# ╔═╡ 756bf2ee-55c1-11eb-2363-d51f3f052b4f
md"""
Effect sizes can be classified as Medium, Insignificant, Small, and Large, respectively.
"""

# ╔═╡ b6c33e70-55c1-11eb-3892-d321937a86ed
md"""
### Check presence of interactions
"""

# ╔═╡ e3220a00-55c1-11eb-080c-451cd462a4d6
plot(medianresult)

# ╔═╡ cbfc9a20-55d0-11eb-3cd3-23a3413a033c
pwd()

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
# ╠═f0393540-5514-11eb-2f76-3514cc3bc85c
# ╟─f02df530-5515-11eb-0d9e-8dbb63ed2023
# ╟─1bf98170-5516-11eb-0c24-77ac43672fba
# ╟─25253910-5516-11eb-0a7a-1d73ed51d4a1
# ╟─d8c96c5e-5517-11eb-0861-056d5e425f04
# ╟─87eb8f70-5518-11eb-1639-f91b4bcb5c9b
# ╟─a87aac80-5518-11eb-38d4-5fe90140b820
# ╟─c8f83d60-5518-11eb-2643-a30df2dd9fb9
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
# ╟─5f1ce3b0-55c1-11eb-3d03-b53516cfc976
# ╟─756bf2ee-55c1-11eb-2363-d51f3f052b4f
# ╟─b6c33e70-55c1-11eb-3892-d321937a86ed
# ╠═e3220a00-55c1-11eb-080c-451cd462a4d6
# ╠═cbfc9a20-55d0-11eb-3cd3-23a3413a033c
