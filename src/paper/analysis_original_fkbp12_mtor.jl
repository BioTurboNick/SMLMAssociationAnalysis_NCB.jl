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
# Analysis of Experiment 3: FKBP12-mTOR
"""

# ╔═╡ 90d2acb0-5507-11eb-39c4-df93c258dd98
md"""Load packages"""

# ╔═╡ 86502ff0-550d-11eb-1c13-553251c061ee
md"""Load results file and compute normalized values"""

# ╔═╡ 5e3db410-550d-11eb-0331-49a447464b07
begin
	### Load saved data

	nreplicates = 5
	nsamples = 2
	ncells = 10

	outputdir = "../../output"
	datapath = joinpath(outputdir, "resultspos.jld2")

	experimentresults = load(datapath)["experimentresults"]

	medianmeasurements = Array{Float64,3}(undef, ncells, nreplicates, nsamples)
	montecarlomeasurements = Array{Float64,3}(undef, ncells, nreplicates, nsamples)
	positivecontrolmontecarlomeasurements = Array{Float64,3}(undef, ncells, nreplicates, nsamples)
	negativecontrolmontecarlomeasurements = Array{Float64,3}(undef, ncells, nreplicates, nsamples)

	for i ∈ 1:nsamples
		samplemedianresults = Array{Float64,2}(undef, ncells, nreplicates)
		samplelessthan10results = Array{Float64,2}(undef, ncells, nreplicates)
		samplepositivecontrollessthan10results = Array{Float64,2}(undef, ncells, nreplicates)
		samplenegativecontrollessthan10results = Array{Float64,2}(undef, ncells, nreplicates)

		for j ∈ 1:nreplicates
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
- 3 - condition 1 (1 = -Rap, 2 = +Rap)
"""

# ╔═╡ 9a94f370-5511-11eb-1f60-7fe83ed11def
md"""
### Planning
Central Limit Theorem - Only 10 measurements per level, but 30 per condition. Might be close enough for a normal distribution.

Group sizes equal
"""

# ╔═╡ 572ba5f2-5513-11eb-1937-3f03cec8252e
md"""
### Medians
#### Check unusual cases
"""

# ╔═╡ 864b181e-5513-11eb-0304-8983d37a8170
let
	p1 = boxplot(medianmeasurements[:, :, 1], xaxis = ("Replicates", [1, 2, 3, 4, 5]), yaxis = ("Median distance (nm)"))
	p2 = boxplot(medianmeasurements[:, :, 2], xaxis = ("Replicates", [1, 2, 3, 4, 5]), yaxis = ("Median distance (nm)"))

	plot(p1, p2, layout = grid(1, 2), legend = :none, plot_title = "FKBP12-mTOR")
end

# ╔═╡ 27e24820-5514-11eb-1415-8711a194892d
md"""
Potential outliers in box plots:
"""

# ╔═╡ f0393540-5514-11eb-2f76-3514cc3bc85c
[medianmeasurements[4,1,1],
 medianmeasurements[[4,9],5,1],
 medianmeasurements[4,2,2],
 medianmeasurements[9,3,2]]

# ╔═╡ f02df530-5515-11eb-0d9e-8dbb63ed2023
md"""Evaluating z-scores within each replicate:"""

# ╔═╡ ffa9f6c0-576e-11eb-0f11-3139a34ff859
mediansflat = [medianmeasurements[:, :, 1] medianmeasurements[:, :, 2]]

# ╔═╡ 25253910-5516-11eb-0a7a-1d73ed51d4a1
medianzscores = [zscore(mediansflat[:,i]) for i ∈ axes(mediansflat, 2)]

# ╔═╡ d8c96c5e-5517-11eb-0861-056d5e425f04
md"""
- In [1] Replicate 1: 1 value, Cell 4 (10%) outside 2.58. None outside 3.29.
- In [1] Replicate 2: All within 1.96
- In [1] Replicate 3: 1 value, Cell 2 (10%) outside 1.96 significantly. None outside 2.58.
- In [1] Replicate 4: All within 1.96
- In [1] Replicate 5: 1 value, Cell 9 (10%) outside 1.96 significantly. None outside 2.58.
- In [2] Replicate 1: All within 1.96
- In [2] Replicate 2: 1 value, Cell 4 (10%) just under 2.58. None outside 2.58.
- In [2] Replicate 3: 1 value, Cell 9 (10%) over 1.96. None outside 2.58.
- In [2] Replicate 4: All within 1.96
- In [2] Replicate 5: All within 1.96
"""

# ╔═╡ 87eb8f70-5518-11eb-1639-f91b4bcb5c9b
md"""
None of these potential outliers are due to data entry error or measurement error. Sampling a different population is possible but unknown. These potential outliers are relatively modest, though one quite severe. But since an ANOVA will be used, must be careful not to delete points from a balanced design. Therefore, no outliers will be removed.

However, those two outliers in [1] Replicate 1, 3, and 5, and [2] Replicate 2 are so extreme that I will windsorize them.
"""

# ╔═╡ 6da539c0-5907-11eb-0f68-c7a5f1ee90ca
begin
	medianmeasurementsw = copy(medianmeasurements)
	mediansflatw = copy(mediansflat)
	medianmeasurementsw[4,1,1] = medianmeasurementsw[5,1,1]
	mediansflatw[4,1] = mediansflatw[5,1]
	medianmeasurementsw[2,3,1] = medianmeasurementsw[5,3,1]
	mediansflatw[2,3] = mediansflatw[5,3]
	medianmeasurementsw[9,5,1] = medianmeasurementsw[4,5,1]
	mediansflatw[9,5] = mediansflatw[4,5]
	medianmeasurementsw[4,2,2] = medianmeasurementsw[10,2,2]
	mediansflatw[4,7] = mediansflatw[10,7]
end

# ╔═╡ a87aac80-5518-11eb-38d4-5fe90140b820
md"""
#### Check assumption of homoscedasticity
"""

# ╔═╡ f6e23cd0-5518-11eb-226a-c5751eb65b34
md"""
Based on the ZResid/ZPred plot, the data appears fairly heteroscedastic; outliers really stand out.

Conducted Levene's test, a 1-way ANOVA of the absolute deviation between each value and the group mean, by replicate.
"""

# ╔═╡ 0ffdd990-5519-11eb-376b-33be0cc966d7
levene(mediansflatw)

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
medianskewness = [skewness(mediansflatw[:,i]) for i ∈ axes(mediansflatw, 2)]

# ╔═╡ 135dbe6e-55be-11eb-3a89-dd3d02c3ad4a
mediankurtosis = [kurtosis(mediansflatw[:,i]) for i ∈ axes(mediansflatw, 2)]

# ╔═╡ 6778a880-55be-11eb-1016-c59335beef0a
md"""
From skewness and kurtosis statistics, all values were within 1.96 (5% level).
"""

# ╔═╡ 9bf9aaf0-55be-11eb-17a2-75b3eec550f8
md"""
The primary source of non-normality appear to be outliers, relatively minimal. While I will further windsorize the outlier in [2] Replicate 2, I won't intend to do anything special.
"""

# ╔═╡ 265006c0-55c1-11eb-2ba8-3d628d096704
md"""
#### Conduct ANOVA tests
"""

# ╔═╡ c101d900-57ab-11eb-04ab-77cd4e0ff4b3
medianresult = anova(
    permutedims(medianmeasurementsw, (1,3,2)),
    [fixed, subject],
    factornames = ["Rapamycin", "Replicate"],
)

# ╔═╡ fd2d4a40-57ab-11eb-2688-53b009a37680
md"""
#### Interpret main effects
Rapamycin does not cause a significant effect on median distance.
"""

# ╔═╡ 895181b0-59f2-11eb-3d4e-1d92833b0e72
md"""
#### Figure
"""

# ╔═╡ 95fd3850-59f2-11eb-0914-a999ed3b9e5f
let
	import StatsPlots.mm
	fkbp12_mtor_median = [medianmeasurementsw[:,:,1] |> vec; medianmeasurementsw[:,:,2] |> vec]
	groups = repeat([1,2], inner = 50)
	boxplot(groups, fkbp12_mtor_median, outliers=false,
			legend=:none,
			bottom_margin=5mm,
			seriescolor=[:white :lightgray],
			line=(6, 0.75),
			xaxis=("Rapamycin", (1:2, ["-Rap", "+Rap"])),
			yaxis=("Median distance (nm)"))
	fkbp12_mtor_median_means = dropdims(mean(medianmeasurementsw, dims=1), dims=1)
	dotplot!([1,1,1,1,1,2,2,2,2,2], fkbp12_mtor_median_means |> vec, mode = :none, label="", marker=(12, 0.75, :rect, repeat([:orange, :darkblue, :darkred, :darkgreen, :darkgray]), stroke(0)))
	dotplot!(groups, fkbp12_mtor_median, mode = :density, label="", marker=(8, 0.5, repeat([:orange, :darkblue, :darkred, :darkgreen, :darkgray], inner=10), stroke(0)))
end

# ╔═╡ 25d0ce1e-5691-11eb-2cc8-49980b6b2161
md"""
### Monte Carlo fraction associated
#### Check unusual cases
"""

# ╔═╡ 8c0ebd10-56a9-11eb-3e8c-95b53bf12061
let
	p1 = boxplot(montecarlomeasurements[:, :, 1], xaxis = ("Replicates", [1, 2, 3, 4, 5]), yaxis = ("Fraction bound"))
	p2 = boxplot(montecarlomeasurements[:, :, 2], xaxis = ("Replicates", [1, 2, 3, 4, 5]), yaxis = ("Fraction bound"))

	plot(p1, p2, layout = grid(1, 2), legend = :none, plot_title = "FKBP12-mTOR")
end

# ╔═╡ faa812c0-56aa-11eb-3b6b-99fe5f200b67
md"""
Potential outliers in box plots:
"""

# ╔═╡ 012a7c50-56ab-11eb-3d50-c9820925f9ac
[montecarlomeasurements[8,[3,4],1],
 montecarlomeasurements[[4,10],1,2],
 montecarlomeasurements[[4,9],2,2],
 montecarlomeasurements[8,5,2]]

# ╔═╡ 1637a040-56ac-11eb-2433-49c37383e005
md"""Evaluating z-scores within each replicate:"""

# ╔═╡ 87ec3520-56ac-11eb-2a1a-2981155d4168
montecarloflat = [montecarlomeasurements[:, :, 1] montecarlomeasurements[:, :, 2]]

# ╔═╡ 9ce1f5f0-56ac-11eb-22a5-d9b5d6a18978
montecarlozscores = [zscore(montecarloflat[:,i]) for i ∈ axes(montecarloflat, 2)]

# ╔═╡ ba7d2760-56ac-11eb-11ca-d543a628c68b
md"""
- In [1] Replicate 1: 1 value, Cell 4 (10%) outside 1.96 barely. None outside 2.58.
- In [1] Replicate 2: All within 1.96
- In [1] Replicate 3: 1 value, Cell 8 (10%) outside 1.96 significantly. None outside 2.58.
- In [1] Replicate 4: 1 value, Cell 8 (10%) outside 1.96 significantly. None outside 2.58.
- In [1] Replicate 5: All within 1.96
- In [2] Replicate 1: All within 1.96
- In [2] Replicate 2: All within 1.96
- In [1] Replicate 3: All within 1.96
- In [2] Replicate 4: All within 1.96
- In [2] Replicate 5: 1 value, Cell 8 (10%) just under 2.58. None outside 2.58.
"""

# ╔═╡ c18cc050-56ad-11eb-261c-4177aa31ca7c
md"""
None of these potential outliers are due to data entry error or measurement error. Sampling a different population is possible but unknown. These potential outliers are relatively modest. In addition, since an ANOVA will be used, must be careful not to delete points from a balanced design. Therefore, no outliers will be removed.

However, I will Windsorize the [1] Replicate 3 and 4, and [2] Replicate 5 outliers.
"""

# ╔═╡ 0498eae0-59e7-11eb-028c-47a1e64c4e5e
begin
	montecarlomeasurementsw = copy(montecarlomeasurements)
	montecarloflatw = copy(montecarloflat)
	montecarlomeasurementsw[8,3,1] = montecarlomeasurementsw[7,3,1]
	montecarlomeasurementsw[8,4,1] = montecarlomeasurementsw[6,4,1]
	montecarlomeasurementsw[8,5,2] = montecarlomeasurementsw[4,5,2]
	montecarloflatw[8,3] = montecarloflatw[7,3]
	montecarloflatw[8,4] = montecarloflatw[6,4]
	montecarloflatw[8,10] = montecarloflatw[4,10]
end

# ╔═╡ 2bba9b02-56ae-11eb-2e30-439dc20273f8
md"""
#### Check assumption of homoscedasticity
"""

# ╔═╡ 0f8a00a0-56af-11eb-2dc6-abd2a45298a4
md"""
Based on the ZResid/ZPred plot, the data appears somewhat heteroscedastic.

Conducted Levene's test, a 1-way ANOVA of the absolute deviation between each value and the group mean, by replicate.
"""

# ╔═╡ b1ae0750-56af-11eb-0714-35b78e197b5e
levene(montecarloflatw)

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
montecarloskewness = [skewness(montecarloflatw[:,i]) for i ∈ axes(montecarloflatw, 2)]

# ╔═╡ 047df370-56b3-11eb-3036-2393031bd967
montecarlokurtosis = [kurtosis(montecarloflatw[:,i]) for i ∈ axes(montecarloflatw, 2)]

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
    permutedims(montecarlomeasurementsw, (1,3,2)),
	[fixed, subject],
    factornames = ["Rapamycin", "Replicate"],
)

# ╔═╡ 2553a3a0-57a9-11eb-247d-a161d85d4410
md"""
#### Interpret main effects
Rapamycin has a large effect and is a bit above the significance threshold.
"""

# ╔═╡ 38571c12-59f3-11eb-337f-13164dbd552b
md"""
#### Figure
"""

# ╔═╡ 475acc70-59f3-11eb-08a1-99e272fbc933
let
	import StatsPlots.mm
	fkbp12_mtor_montecarlo = [montecarlomeasurementsw[:,:,1] |> vec; montecarlomeasurementsw[:,:,2] |> vec]
	groups = repeat([1,2], inner = 50)
	boxplot(groups, fkbp12_mtor_montecarlo, outliers=false,
			legend=:none,
			bottom_margin=5mm,
			seriescolor=[:white :lightgray],
			line=(6, 0.75),
		size=(512, 1024),
			xaxis=("Rapamycin", (1:2, ["-Rap", "+Rap"])),
			yaxis=("Fraction associated"))
	fkbp12_mtor_montecarlo_means = dropdims(mean(montecarlomeasurementsw, dims=1), dims=1)
	dotplot!([1,1,1,1,1,2,2,2,2,2], fkbp12_mtor_montecarlo_means |> vec, mode = :none, label="", marker=(12, 0.75, :rect, repeat([:orange, :darkblue, :darkred, :darkgreen, :darkgray]), stroke(0)))
	dotplot!(groups, fkbp12_mtor_montecarlo, mode = :density, label="", marker=(8, 0.5, repeat([:orange, :darkblue, :darkred, :darkgreen, :darkgray], inner=10), stroke(0)))
end

# ╔═╡ 78883b10-5756-11eb-32cf-b95fdeebda79
md"""
### Normalized Monte Carlo fraction associated
This approach was abandoned because we determined it was over-normalizing the data, but still included here for completeness.

#### Check unusual cases
"""

# ╔═╡ df7dce20-5756-11eb-3796-037e8e3c49f9
let
	p1 = boxplot(normalizedmontecarlomeasurements[:, :, 1], xaxis = ("Replicates", [1, 2, 3, 4, 5]), yaxis = ("Fraction bound"))
	p2 = boxplot(normalizedmontecarlomeasurements[:, :, 2], xaxis = ("Replicates", [1, 2, 3, 4, 5]), yaxis = ("Fraction bound"))

	plot(p1, p2, layout = grid(1, 2), legend = :none, plot_title = "FKBP12-mTOR")
end

# ╔═╡ faf5aa60-5756-11eb-3b5b-f1ad5e004af2
md"""
Potential outliers in box plots:
"""

# ╔═╡ ff845a40-5756-11eb-2ce8-d5ba8458e2ae
[normalizedmontecarlomeasurements[[1,8],4,1],
 normalizedmontecarlomeasurements[8,3,1],
 normalizedmontecarlomeasurements[[3,5,8],4,1],
 normalizedmontecarlomeasurements[6,5,1],
 normalizedmontecarlomeasurements[[4,10],5,1],]

# ╔═╡ d2282670-5757-11eb-0bd9-b1aaf0b37d70
md"""Evaluating z-scores within each replicate:"""

# ╔═╡ d99fec80-5757-11eb-39cd-4b0febbf61e3
normmontecarloflat = [normalizedmontecarlomeasurements[:, :, 1] normalizedmontecarlomeasurements[:, :, 2]]

# ╔═╡ edee22b2-5757-11eb-317a-d3e4d31486cd
normmontecarlozscores = [zscore(normmontecarloflat[:,i]) for i ∈ axes(normmontecarloflat, 2)]

# ╔═╡ a05efb90-5758-11eb-10ca-fb9fb8e05821
md"""
- In [1] Replicate 1: All within 1.96
- In [1] Replicate 2: All within 1.96
- In [1] Replicate 3: 1 value, Cell 8 (10%) outside 1.96. None outside 2.58.
- In [1] Replicate 4: 1 value, Cell 8 (10%) outside 1.96 significantly. None outside 2.58.
- In [1] Replicate 5: 1 value, Cell 6 (10%) outside 1.96 significantly. None outside 2.58.
- In [2] Replicate 1: All within 1.96
- In [2] Replicate 2: All within 1.96
- In [2] Replicate 3: All within 1.96
- In [2] Replicate 4: All within 1.96
- In [2] Replicate 5: 1 value, Cell 8 (10%) outside 1.96. None outside 2.58.
"""

# ╔═╡ 7323e450-5759-11eb-1c62-7ba711bf418b
md"""
None of these potential outliers are due to data entry error or measurement error. Sampling a different population is possible but unknown. These potential outliers are relatively modest. In addition, since an ANOVA will be used, must be careful not to delete points from a balanced design. Therefore, no outliers will be removed.

However, I will Windsorize the worst outliers of [1] Replicate 4 and 5.
"""

# ╔═╡ 49094280-66a6-11eb-1261-7fe254f02b6d
begin
	normalizedmontecarlomeasurementsw = copy(normalizedmontecarlomeasurements)
	normmontecarloflatw = copy(normmontecarloflat)
	normalizedmontecarlomeasurementsw[8,4,1] = normalizedmontecarlomeasurementsw[5,4,1]
	normalizedmontecarlomeasurementsw[6,5,1] = normalizedmontecarlomeasurementsw[9,5,1]
	normmontecarloflatw[8,4] = normmontecarloflatw[5,4]
	normmontecarloflatw[6,5] = normmontecarloflatw[9,5]
end

# ╔═╡ df4eee90-5759-11eb-02ef-53a0c6f31f07
md"""
#### Check assumption of homoscedasticity
"""

# ╔═╡ 8c2a8790-575b-11eb-34be-99f584826d43
md"""
Based on the ZResid/ZPred plot, the data appears somewhat heteroscedastic, although not much relationship between residual and predicted.

Conducted Levene's test, a 1-way ANOVA of the absolute deviation between each value and the group mean, by replicate.
"""

# ╔═╡ 9aade460-575b-11eb-29da-4ff6881ceb60
levene(normmontecarloflatw)

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
normmontecarloskewness = [skewness(normmontecarloflatw[:,i]) for i ∈ axes(normmontecarloflatw, 2)]

# ╔═╡ 682c7090-575d-11eb-3b12-d50d450d3fce
normmontecarlokurtosis = [kurtosis(normmontecarloflatw[:,i]) for i ∈ axes(normmontecarloflatw, 2)]

# ╔═╡ 80f22522-575d-11eb-211b-59fe55838280
md"""
From skewness and kurtosis statistics, one was greater in value than 1.96 (5% level), [1] Replicate 4, with skewness 1.78 and kurtosis 2.36.
"""

# ╔═╡ 91c676d0-575d-11eb-2fd5-97ca14bfdec0
md"""
#### Conduct ANOVA tests
"""

# ╔═╡ b7624d5e-575d-11eb-1221-974537c3b48e
normmontecarloresult = anova(
    permutedims(normalizedmontecarlomeasurementsw, (1, 3, 2)),
    [fixed, subject],
    factornames = ["Rapamycin", "Replicate"],
)

# ╔═╡ 408e5d30-5764-11eb-124f-4f2630d1d908
md"""
#### Interpret main effects
Rapamycin had a moderate effect but was not significant (2.48% to 5.85%).
"""

# ╔═╡ 8c39ec6e-5a72-11eb-3b9b-1d45fbdeaa8c
md"""
#### Figure
"""

# ╔═╡ 9cd79dc0-5a72-11eb-12a5-cd1e42c92297
let
	import StatsPlots.mm
	fkbp12_mtor_montecarlo = [normalizedmontecarlomeasurements[:,:,1] |> vec; normalizedmontecarlomeasurements[:,:,2] |> vec]
	groups = repeat([1,2], inner = 50)
	boxplot(groups, fkbp12_mtor_montecarlo, outliers=false,
			legend=:none,
			bottom_margin=5mm,
			seriescolor=[:white :lightgray],
			line=(6, 0.75),
			xaxis=("Rapamycin", (1:2, ["-Rap", "+Rap"])),
			yaxis=("Median distance (nm)"))
	fkbp12_mtor_montecarlo_means = dropdims(mean(normalizedmontecarlomeasurements, dims=1), dims=1)
	dotplot!([1,1,1,1,1,2,2,2,2,2], fkbp12_mtor_montecarlo_means |> vec, mode = :none, label="", marker=(12, 0.75, :rect, repeat([:orange, :darkblue, :darkred, :darkgreen, :darkgray]), stroke(0)))
	dotplot!(groups, fkbp12_mtor_montecarlo, mode = :density, label="", marker=(8, 0.5, repeat([:orange, :darkblue, :darkred, :darkgreen, :darkgray], inner=10), stroke(0)))
end

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
zresid_zpred_plot(mediansflatw)

# ╔═╡ 3a3fa580-56ae-11eb-1f19-f3c06ef2d9d7
zresid_zpred_plot(montecarloflatw)

# ╔═╡ 59522180-575a-11eb-39c4-f9eaa27f7901
zresid_zpred_plot(normmontecarloflatw)

# ╔═╡ 7922d1d0-56b0-11eb-261c-5f87ca587394
function qqnormplot(data)
	zresid, _ = zresid_zpred(data)
	p = [qqnorm(zresid[:, i]) for i ∈ axes(data, 2)]
	plot(p..., layout = grid(5, 2), legend = :none, plot_title = "MDM2-p53 qqnorm")
end

# ╔═╡ 6f727ed0-5519-11eb-1462-01c83aa975b6
qqnormplot(mediansflatw)

# ╔═╡ 5ef3f050-56b0-11eb-31f9-419a4161b67e
qqnormplot(montecarloflatw)

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
# ╠═f0393540-5514-11eb-2f76-3514cc3bc85c
# ╟─f02df530-5515-11eb-0d9e-8dbb63ed2023
# ╟─ffa9f6c0-576e-11eb-0f11-3139a34ff859
# ╟─25253910-5516-11eb-0a7a-1d73ed51d4a1
# ╟─d8c96c5e-5517-11eb-0861-056d5e425f04
# ╟─87eb8f70-5518-11eb-1639-f91b4bcb5c9b
# ╠═6da539c0-5907-11eb-0f68-c7a5f1ee90ca
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
# ╠═c101d900-57ab-11eb-04ab-77cd4e0ff4b3
# ╟─fd2d4a40-57ab-11eb-2688-53b009a37680
# ╟─895181b0-59f2-11eb-3d4e-1d92833b0e72
# ╟─95fd3850-59f2-11eb-0914-a999ed3b9e5f
# ╟─25d0ce1e-5691-11eb-2cc8-49980b6b2161
# ╟─8c0ebd10-56a9-11eb-3e8c-95b53bf12061
# ╟─faa812c0-56aa-11eb-3b6b-99fe5f200b67
# ╠═012a7c50-56ab-11eb-3d50-c9820925f9ac
# ╟─1637a040-56ac-11eb-2433-49c37383e005
# ╟─87ec3520-56ac-11eb-2a1a-2981155d4168
# ╟─9ce1f5f0-56ac-11eb-22a5-d9b5d6a18978
# ╟─ba7d2760-56ac-11eb-11ca-d543a628c68b
# ╟─c18cc050-56ad-11eb-261c-4177aa31ca7c
# ╠═0498eae0-59e7-11eb-028c-47a1e64c4e5e
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
# ╟─38571c12-59f3-11eb-337f-13164dbd552b
# ╟─475acc70-59f3-11eb-08a1-99e272fbc933
# ╟─78883b10-5756-11eb-32cf-b95fdeebda79
# ╟─df7dce20-5756-11eb-3796-037e8e3c49f9
# ╟─faf5aa60-5756-11eb-3b5b-f1ad5e004af2
# ╠═ff845a40-5756-11eb-2ce8-d5ba8458e2ae
# ╟─d2282670-5757-11eb-0bd9-b1aaf0b37d70
# ╟─d99fec80-5757-11eb-39cd-4b0febbf61e3
# ╠═edee22b2-5757-11eb-317a-d3e4d31486cd
# ╟─a05efb90-5758-11eb-10ca-fb9fb8e05821
# ╟─7323e450-5759-11eb-1c62-7ba711bf418b
# ╠═49094280-66a6-11eb-1261-7fe254f02b6d
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
# ╟─408e5d30-5764-11eb-124f-4f2630d1d908
# ╟─8c39ec6e-5a72-11eb-3b9b-1d45fbdeaa8c
# ╟─9cd79dc0-5a72-11eb-12a5-cd1e42c92297
# ╟─9c92b290-56ae-11eb-2595-85845fe95f0e
# ╟─2623d480-56af-11eb-2ab8-0d0afd8c45dc
# ╟─a655a120-56ae-11eb-30be-459d7c3067dc
# ╟─7922d1d0-56b0-11eb-261c-5f87ca587394
