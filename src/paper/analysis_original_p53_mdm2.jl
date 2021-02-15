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
# Analysis of Experiment 1: p53-Mdm2
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

	nothing
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

Doxycycline and Nutlin treatments should be independent.
"""

# ╔═╡ 572ba5f2-5513-11eb-1937-3f03cec8252e
md"""
### Medians
#### Check unusual cases
"""

# ╔═╡ 864b181e-5513-11eb-0304-8983d37a8170
let
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
- In [2,1] Replicate 3: 1 value, Cell 10 (10%) outside 1.96 significantly. None outside 2.58.
- In [2,2] Replicate 1: All within 1.96
- In [2,2] Replicate 2: All within 1.96
- In [2,2] Replicate 3: 1 value, Cell 5 (10%) outside 1.96 significantly. None outside 2.58.
"""

# ╔═╡ 87eb8f70-5518-11eb-1639-f91b4bcb5c9b
md"""
None of these potential outliers are due to data entry error or measurement error. Sampling a different population is possible but unknown. These potential outliers are relatively modest. In addition, since an ANOVA will be used, must be careful not to delete points from a balanced design. Therefore, no outliers will be removed. However, I will Windsorize the more extreme ones for [1,2] Replicate 2, [2,1] Replicate 3, and [2,2] Replicate 3.
"""

# ╔═╡ 9b7f65a2-5952-11eb-0c21-733df3703bf9
begin
	medianmeasurementsw = copy(medianmeasurements)
	mediansflatw = copy(mediansflat)
	medianmeasurementsw[6,2,2,1,1] = medianmeasurements[7,2,2,1,1]
	medianmeasurementsw[10,3,1,2,1] = medianmeasurements[4,3,1,2,1]
	medianmeasurementsw[5,3,2,2,1] = medianmeasurements[9,3,2,2,1]
	mediansflatw[6,5] = mediansflat[7,5]
	mediansflatw[10,9] = mediansflat[4,9]
	mediansflatw[5,12] = mediansflat[9,12]
end;

# ╔═╡ a87aac80-5518-11eb-38d4-5fe90140b820
md"""
#### Check assumption of homoscedasticity
"""

# ╔═╡ f6e23cd0-5518-11eb-226a-c5751eb65b34
md"""
Based on the ZResid/ZPred plot, the data appears fairly heteroscedastic.

Conducted Levene's test, a 1-way ANOVA of the absolute deviation between each value and the group mean, by replicate.
"""

# ╔═╡ 0ffdd990-5519-11eb-376b-33be0cc966d7
levene(mediansflatw)

# ╔═╡ 28a8da30-5519-11eb-2fb0-efe693c6e44c
md"""
By Levene's test, the data is very likely heteroscedastic. Based on median, still significant at p = 0.011.
 
Not doing Hartley's Fmax because it isn't as simple.
 
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

# ╔═╡ 5f1ce3b0-55c1-11eb-3d03-b53516cfc976
medianresult = anova(
    permutedims(medianmeasurementsw[:, :, :, :, 1], (1,3,4,2)),
    [fixed, fixed, subject],
    factornames = ["Doxycycline", "Nutlin-3a", "Replicate"],
)

# ╔═╡ 756bf2ee-55c1-11eb-2363-d51f3f052b4f
md"""
Effect sizes can be classified as Medium, Insignificant, Small, and Large, respectively.
"""

# ╔═╡ b6c33e70-55c1-11eb-3892-d321937a86ed
md"""
#### Check presence of interactions
"""

# ╔═╡ e3220a00-55c1-11eb-080c-451cd462a4d6
plot(medianresult)

# ╔═╡ cbfc9a20-55d0-11eb-3cd3-23a3413a033c
md"""
The lines are quite different, suggesting a potential strong interaction. However, the F statistic for the interaction was not significant. Perhaps the number of replicates isn't enough due to variance between replicates. May be worth checking simple effects anyway.
"""

# ╔═╡ 2ed8ff90-568e-11eb-10a3-0b97b2425d2f
md"""
#### Interpret main effects
Nutlin-3a is just above the significance threshold with a medium effect size while Doxycyline produced no changes whatsoever. Nutlin-3a was thus associated with a drop in median distance, likely due to an increase in expression level.
"""

# ╔═╡ 76abe9e0-6af2-11eb-0dc5-4594d89c209b
md"""
### Figure
"""

# ╔═╡ 81b932c0-6af2-11eb-0a5c-09c0b3f10519
let
	import StatsPlots.mm
	p53_mdm2_median = [medianmeasurementsw[:,:,1,1,1] |> vec; medianmeasurementsw[:,:,2,1,1] |> vec; medianmeasurementsw[:,:,1,2,1] |> vec; medianmeasurementsw[:,:,2,2,1] |> vec]
	groups = repeat([2, 1, 4, 3], inner = 30)
	boxplot(groups, p53_mdm2_median, outliers=false,
		    label=["- Dox", "+ Dox"],
			legend=:none,
			seriescolor=:white,
			size=(256, 512),
			line=(2, 0.75),
		    xgrid=:none,
			xaxis=("Treatment", (1:4, ["-Dox\n-Nut", "+Dox\n-Nut", "-Dox\n+Nut", "+Dox\n+Nut"])),
			yaxis=("Median distance (nm)", [0, 600]))
	p53_mdm2_median_means = dropdims(mean(medianmeasurementsw[:,:,:,:,1], dims = 1), dims = 1) |> vec
	dotplot!(repeat([2, 1, 4, 3], inner=3), p53_mdm2_median_means |> vec, mode = :none, label="", marker=(4, 0.75, :rect, repeat([:orange, :darkblue, :darkred]), stroke(0)))
	dotplot!(groups, p53_mdm2_median, mode = :density, label="", marker=(2, 0.5, repeat([:orange, :darkblue, :darkred], inner=10), stroke(0)))
end

# ╔═╡ b69b87b2-6af5-11eb-30c2-bfa8159a8e0d
md"""Create paper-quality version and save it:"""

# ╔═╡ deafcdc0-6af4-11eb-1673-0101d277b4b3
let
	import StatsPlots.mm
	p53_mdm2_median = [medianmeasurementsw[:,:,1,1,1] |> vec; medianmeasurementsw[:,:,2,1,1] |> vec; medianmeasurementsw[:,:,1,2,1] |> vec; medianmeasurementsw[:,:,2,2,1] |> vec]
	groups = repeat([2, 1, 4, 3], inner = 30)
	boxplot(groups, p53_mdm2_median, outliers=false,
			legend=:none,
		    left_margin=25mm,
			seriescolor=:white,
		    guidefontsize=52,
            tickfontsize=48,
			size=(1024, 2048),
			line=(8, 0.75),
		    xgrid=:none,
			xaxis=("Treatment", (1:4, ["-Dox\n-Nut", "+Dox\n-Nut", "-Dox\n+Nut", "+Dox\n+Nut"])),
			yaxis=("Median distance (nm)", [0, 600]))
	p53_mdm2_median_means = dropdims(mean(medianmeasurementsw[:,:,:,:,1], dims=1), dims=1) |> vec
	dotplot!(repeat([2, 1, 4, 3], inner=3), p53_mdm2_median_means |> vec, mode = :none, label="", marker=(12, 0.75, :rect, repeat([:orange, :darkblue, :darkred]), stroke(0)))
	dotplot!(groups, p53_mdm2_median, mode = :density, label="", marker=(8, 0.5, repeat([:orange, :darkblue, :darkred], inner=10), stroke(0)))
	savefig(joinpath(outputdir, "p53-mdm2-medians.png"))
end

# ╔═╡ 25d0ce1e-5691-11eb-2cc8-49980b6b2161
md"""
### Monte Carlo fraction associated
#### Check unusual cases
"""

# ╔═╡ 8c0ebd10-56a9-11eb-3e8c-95b53bf12061
let
	p1 = boxplot(montecarlomeasurements[:, :, 1, 1, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p2 = boxplot(montecarlomeasurements[:, :, 2, 1, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p3 = boxplot(montecarlomeasurements[:, :, 1, 2, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p4 = boxplot(montecarlomeasurements[:, :, 2, 2, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))

	plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "Mdm2-p53")
end

# ╔═╡ faa812c0-56aa-11eb-3b6b-99fe5f200b67
md"""
Potential outliers in box plots:
"""

# ╔═╡ 012a7c50-56ab-11eb-3d50-c9820925f9ac
[montecarlomeasurements[6,1,2,1,1],
 montecarlomeasurements[1,3,1,2,1],
 montecarlomeasurements[5,3,2,2,1]]

# ╔═╡ 1637a040-56ac-11eb-2433-49c37383e005
md"""Evaluating z-scores within each replicate:"""

# ╔═╡ 87ec3520-56ac-11eb-2a1a-2981155d4168
montecarloflat = [montecarlomeasurements[:, :, 1, 1, 1] montecarlomeasurements[:, :, 2, 1, 1] montecarlomeasurements[:, :, 1, 2, 1] montecarlomeasurements[:, :, 2, 2, 1]]

# ╔═╡ 9ce1f5f0-56ac-11eb-22a5-d9b5d6a18978
montecarlozscores = [zscore(montecarloflat[:,i]) for i ∈ axes(montecarloflat, 2)]

# ╔═╡ ba7d2760-56ac-11eb-11ca-d543a628c68b
md"""
- In [1,1] Replicate 1: All within 1.96
- In [1,1] Replicate 2: All within 1.96
- In [1,1] Replicate 3: All within 1.96
- In [1,2] Replicate 1: 1 value, Cell 6 (10%) outside 1.96 significantly. None outside 2.58.
- In [1,2] Replicate 2: All within 1.96
- In [1,2] Replicate 3: All within 1.96
- In [2,1] Replicate 1: All within 1.96
- In [2,1] Replicate 2: All within 1.96
- In [2,1] Replicate 3: All within 1.96
- In [2,2] Replicate 1: 1 value, Cell 7 (10%) outside 1.96. None outside 2.58.
- In [2,2] Replicate 2: 1 value, Cell 9 (10%) outside 1.96. None outside 2.58.
- In [2,2] Replicate 3: All within 1.96
"""

# ╔═╡ c18cc050-56ad-11eb-261c-4177aa31ca7c
md"""
None of these potential outliers are due to data entry error or measurement error. Sampling a different population is possible but unknown. These potential outliers are relatively modest. In addition, since an ANOVA will be used, must be careful not to delete points from a balanced design. Therefore, no outliers will be removed. However, I shall Windsorize the significant one in [1,2] Replicate 1.
"""

# ╔═╡ 935aa530-6d73-11eb-2012-4d5e8b594f42
begin
	montecarlomeasurementsw = copy(montecarlomeasurements)
	montecarloflatw = copy(montecarloflat)
	montecarlomeasurementsw[6,1,2,1,1] = montecarlomeasurementsw[2,1,2,1,1]
	montecarloflatw[6,4] = montecarloflat[2,4]
end;

# ╔═╡ 2bba9b02-56ae-11eb-2e30-439dc20273f8
md"""
#### Check assumption of homoscedasticity
"""

# ╔═╡ 0f8a00a0-56af-11eb-2dc6-abd2a45298a4
md"""
Based on the ZResid/ZPred plot, the data appears fairly heteroscedastic.

Conducted Levene's test, a 1-way ANOVA of the absolute deviation between each value and the group mean, by replicate.
"""

# ╔═╡ b1ae0750-56af-11eb-0714-35b78e197b5e
levene(montecarloflatw)

# ╔═╡ f9676320-56af-11eb-1481-ef9502edded8
md"""
By Levene's test, the data is likely heteroscedastic.

Group sizes are equal, but samples are on the smaller side. CIs and significance tests may be affected.
"""

# ╔═╡ 4e36f730-56b0-11eb-21e8-0d2e7d0a2562
md"""
#### Check assumption of normality
"""

# ╔═╡ 4423b3d0-56b2-11eb-22be-fd86deabef7f
md"""
There are a couple notable deviations but for the *most* part they appear close to normal. [2,2] replicate 1 is the most deviant.
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

# ╔═╡ 755d8510-56b3-11eb-3a0c-53b643ba2053
montecarloresult = anova(
    permutedims(montecarlomeasurementsw[:, :, :, :, 1], (1,3,4,2)),
    [fixed, fixed, subject],
    factornames = ["Doxycycline", "Nutlin-3a", "Replicate"],
)

# ╔═╡ bf915350-574e-11eb-2f66-236c4ec9e600
md"""
Effect sizes can be classified as Large, Small, Insignificant, and Large, respectively.
"""

# ╔═╡ d79abcc0-574e-11eb-3777-4f8b85036436
md"""
#### Check presence of interactions
"""

# ╔═╡ f5f23950-574e-11eb-02b5-2fb9b6272ce6
plot(montecarloresult)

# ╔═╡ 27747350-5752-11eb-12b9-dfffc84b15c8
md"""
Lines are nearly parallel and the interaction F statistic was solidly non-significant. Checking the simple effects is valid.
"""

# ╔═╡ 248d1560-5753-11eb-2156-89cfb6742205
md"""
#### Interpret main effects
Nutlin-3a is below the significance threshold with a large effect size while Doxycyline produced a small effect with low probability. Nutlin-3a was thus associated with a decrease in fraction associated, while MEG3 expression likely was not associated with a change.
"""

# ╔═╡ 552d3550-6af5-11eb-2fcb-fb755145a298
md"""
### Figure
"""

# ╔═╡ 5ab31670-6af5-11eb-2e08-83b5228005d1
let
	import StatsPlots.mm
	p53_mdm2_montecarlo = [montecarlomeasurementsw[:,:,1,1,1] |> vec; montecarlomeasurementsw[:,:,2,1,1] |> vec; montecarlomeasurementsw[:,:,1,2,1] |> vec; montecarlomeasurementsw[:,:,2,2,1] |> vec]
	groups = repeat([2, 1, 4, 3], inner = 30)
	boxplot(groups, p53_mdm2_montecarlo, outliers=false,
		    label=["- Dox", "+ Dox"],
			legend=:none,
			seriescolor=:white,
		    left_margin=2mm,
			size=(256, 512),
			line=(2, 0.75),
		    xgrid=:none,
			xaxis=("Treatment", (1:4, ["-Dox\n-Nut", "+Dox\n-Nut", "-Dox\n+Nut", "+Dox\n+Nut"])),
			yaxis=("Median distance (nm)", [0, 0.25]))
	p53_mdm2_montecarlo_means = dropdims(mean(montecarlomeasurementsw[:,:,:,:,1], dims = 1), dims = 1) |> vec
	dotplot!(repeat([2, 1, 4, 3], inner=3), p53_mdm2_montecarlo_means |> vec, mode = :none, label="", marker=(4, 0.75, :rect, repeat([:orange, :darkblue, :darkred]), stroke(0)))
	dotplot!(groups, p53_mdm2_montecarlo, mode = :density, label="", marker=(2, 0.5, repeat([:orange, :darkblue, :darkred], inner=10), stroke(0)))
end

# ╔═╡ bc4851c0-6af5-11eb-0db5-cb71d67f93e3
md"""Create paper-quality version and save it:"""

# ╔═╡ af04a0e0-6af5-11eb-05ed-99c27fd1a0ef
let
	import StatsPlots.mm
	p53_mdm2_montecarlo = [montecarlomeasurementsw[:,:,1,1,1] |> vec; montecarlomeasurementsw[:,:,2,1,1] |> vec; montecarlomeasurementsw[:,:,1,2,1] |> vec; montecarlomeasurementsw[:,:,2,2,1] |> vec]
	groups = repeat([2, 1, 4, 3], inner = 30)
	boxplot(groups, p53_mdm2_montecarlo, outliers=false,
		    label=["- Dox", "+ Dox"],
			legend=:none,
			seriescolor=:white,
		    left_margin=25mm,
		    guidefontsize=52,
            tickfontsize=48,
			size=(1024, 2048),
			line=(8, 0.75),
		    xgrid=:none,
			xaxis=("Treatment", (1:4, ["-Dox\n-Nut", "+Dox\n-Nut", "-Dox\n+Nut", "+Dox\n+Nut"])),
			yaxis=("Fraction associated", [0, 0.25]))
	p53_mdm2_montecarlo_means = dropdims(mean(montecarlomeasurementsw[:,:,:,:,1], dims=1), dims=1) |> vec
	dotplot!(repeat([2, 1, 4, 3], inner=3), p53_mdm2_montecarlo_means |> vec, mode = :none, label="", marker=(12, 0.75, :rect, repeat([:orange, :darkblue, :darkred]), stroke(0)))
	dotplot!(groups, p53_mdm2_montecarlo, mode = :density, label="", marker=(8, 0.5, repeat([:orange, :darkblue, :darkred], inner=10), stroke(0)))
	savefig(joinpath(outputdir, "p53-mdm2-montecarlo.png"))
end

# ╔═╡ 78883b10-5756-11eb-32cf-b95fdeebda79
md"""
### Normalized Monte Carlo fraction associated
This approach was abandoned because we determined it was over-normalizing the data, but still included here for completeness.

#### Check unusual cases
"""

# ╔═╡ df7dce20-5756-11eb-3796-037e8e3c49f9
let
	p1 = boxplot(normalizedmontecarlomeasurements[:, :, 1, 1, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p2 = boxplot(normalizedmontecarlomeasurements[:, :, 2, 1, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p3 = boxplot(normalizedmontecarlomeasurements[:, :, 1, 2, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))
	p4 = boxplot(normalizedmontecarlomeasurements[:, :, 2, 2, 1], xaxis = ("Replicates", [1, 2, 3]), yaxis = ("Fraction bound"))

	plot(p1, p2, p3, p4, layout = grid(2, 2), legend = :none, plot_title = "Mdm2-p53")
end

# ╔═╡ faf5aa60-5756-11eb-3b5b-f1ad5e004af2
md"""
Potential outliers in box plots:
"""

# ╔═╡ ff845a40-5756-11eb-2ce8-d5ba8458e2ae
[normalizedmontecarlomeasurements[6,1,2,1,1],
 normalizedmontecarlomeasurements[[2,9],3,1,2,1]]

# ╔═╡ d2282670-5757-11eb-0bd9-b1aaf0b37d70
md"""Evaluating z-scores within each replicate:"""

# ╔═╡ d99fec80-5757-11eb-39cd-4b0febbf61e3
normmontecarloflat = [normalizedmontecarlomeasurements[:, :, 1, 1, 1] normalizedmontecarlomeasurements[:, :, 2, 1, 1] normalizedmontecarlomeasurements[:, :, 1, 2, 1] normalizedmontecarlomeasurements[:, :, 2, 2, 1]]

# ╔═╡ edee22b2-5757-11eb-317a-d3e4d31486cd
normmontecarlozscores = [zscore(normmontecarloflat[:,i]) for i ∈ axes(normmontecarloflat, 2)]

# ╔═╡ a05efb90-5758-11eb-10ca-fb9fb8e05821
md"""
- In [1,1] Replicate 1: All within 1.96
- In [1,1] Replicate 2: All within 1.96
- In [1,1] Replicate 3: 1 value, Cell 2 (10%) outside 1.96 barely. None outside 2.58
- In [1,2] Replicate 1: 1 value, Cell 6 (10%) outside 1.96 significantly. None outside 2.58.
- In [1,2] Replicate 2: All within 1.96
- In [1,2] Replicate 3: All within 1.96
- In [2,1] Replicate 1: All within 1.96
- In [2,1] Replicate 2: All within 1.96
- In [2,1] Replicate 3: 1 value, Cell 8 (10%) outside 1.96. None outside 2.58.
- In [2,2] Replicate 1: All within 1.96
- In [2,2] Replicate 2: All within 1.96
- In [2,2] Replicate 3: All within 1.96
"""

# ╔═╡ 7323e450-5759-11eb-1c62-7ba711bf418b
md"""
None of these potential outliers are due to data entry error or measurement error. Sampling a different population is possible but unknown. These potential outliers are relatively modest. In addition, since an ANOVA will be used, must be careful not to delete points from a balanced design. Therefore, no outliers will be removed. However, I will Windsorize the extreme outlier in [1,2] Replicate 1.
"""

# ╔═╡ 9ecc38b0-6d74-11eb-2d1a-ebfa12f11a64
begin
	normalizedmontecarlomeasurementsw = copy(normalizedmontecarlomeasurements)
	normmontecarloflatw = copy(normmontecarloflat)
	normalizedmontecarlomeasurementsw[6,1,2,1,1] = normalizedmontecarlomeasurementsw[4,1,2,1,1]
	normmontecarloflatw[6,4] = normmontecarloflat[4,4]
end;

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
"""

# ╔═╡ b7624d5e-575d-11eb-1221-974537c3b48e
normmontecarloresult = anova(
    permutedims(normalizedmontecarlomeasurementsw[:, :, :, :, 1], (1,3,4,2)),
    [fixed, fixed, subject],
    factornames = ["Doxycycline", "Nutlin-3a", "Replicate"],
)

# ╔═╡ 8f0dad80-575f-11eb-2312-adae92d50e83
md"""
### Check presence of interactions
"""

# ╔═╡ a5489880-575f-11eb-1afa-c3bfaf5dccc8
plot(normmontecarloresult)

# ╔═╡ 4b211122-575f-11eb-0aa9-59fc4e9a2df7
md"""
The lines cross, but the values are so small that it's likely nothing. F-statistic is not significant. Interpreting simple effects valid.
"""

# ╔═╡ 408e5d30-5764-11eb-124f-4f2630d1d908
md"""
#### Interpret main effects
Nutlin-3a did not have a significant effect on fraction bound (from 1.7% to 1.7%). This is unexpected for an inhibitor of Mdm2-p53 binding, but it may reflect that p53 and Mdm2 levels both increase until a steady state is reached, maintaining a similar level of association.

Doxycycline did not have a significant effect, though with a small effect size, though may be due to lack of power in the experiment. There was still a consistent increase with doxycycline exposure (from 1.4% to 1.9%). Thus, MEG3 overexpression does not cause significantly less binding between p53 and Mdm2.
"""

# ╔═╡ d83ab120-6be5-11eb-0b44-176f26137172
md"""
## Localization/molecule maps
"""

# ╔═╡ e5020610-6be5-11eb-25ab-adb9479d6b16
let
	# Exp 2 (p53-Mdm2)
	# A8
	insetx, insety = [15100, 16500], [30000, 31400]
	r = experimentresults[1][3][1][8]
	localizationsplot(r, insetbox = [insetx, insety], forprint = true)
	savefig(joinpath(outputdir, "2 - A8 dSTORM.png"))
	insetplot(r, insetx, insety, include_scalebar = false, forprint = true)
	savefig(joinpath(outputdir, "2 - A8 dSTORM 1400nm.png"))
	p1 = localizationsplot(r, insetbox = [insetx, insety], include_scalebar = true)
	p2 = insetplot(r, insetx, insety, include_scalebar = true)
	plot(p1, p2, layout = grid(1,2), size=(1024, 512), fmt = :png)
end

# ╔═╡ 840e9e20-6be7-11eb-2b80-37ffe2f4495d
let
	# Exp 2 (p53-Mdm2)
	# B5
	insetx, insety = [16050, 17450], [21100, 22500]
	r = experimentresults[1][3][2][5]
	localizationsplot(r, insetbox = [insetx, insety], forprint = true)
	savefig(joinpath(outputdir, "2 - B5 dSTORM.png"))
	insetplot(r, insetx, insety, include_scalebar = false, forprint = true)
	savefig(joinpath(outputdir, "2 - B5 dSTORM 1400nm.png"))
	p1 = localizationsplot(r, insetbox = [insetx, insety], include_scalebar = true)
	p2 = insetplot(r, insetx, insety, include_scalebar = true)
	plot(p1, p2, layout = grid(1,2), size=(1024, 512), fmt = :png)
end

# ╔═╡ 7affa0d0-6be8-11eb-3395-555814cab9be
let
	# Exp 2 (p53-Mdm2)
	# C9
	insetx, insety = [16400, 17800], [22697, 24096]
	r = experimentresults[1][2][3][9]
	localizationsplot(r, insetbox = [insetx, insety], forprint = true)
	savefig(joinpath(outputdir, "2 - C9 dSTORM.png"))
	insetplot(r, insetx, insety, include_scalebar = false, forprint = true)
	savefig(joinpath(outputdir, "2 - C9 dSTORM 1400nm.png"))
	p1 = localizationsplot(r, insetbox = [insetx, insety], include_scalebar = true)
	p2 = insetplot(r, insetx, insety, include_scalebar = true)
	plot(p1, p2, layout = grid(1,2), size=(1024, 512), fmt = :png)
end

# ╔═╡ 13ecfdf0-6bea-11eb-3af4-8d627cccdd2b
let
	# Exp 2 (p53-Mdm2)
	# D4
	insetx, insety = [24100, 25500], [14000, 15400]
	r = experimentresults[1][2][4][4]
	localizationsplot(r, insetbox = [insetx, insety], forprint = true)
	savefig(joinpath(outputdir, "2 - D4 dSTORM.png"))
	insetplot(r, insetx, insety, include_scalebar = false, forprint = true)
	savefig(joinpath(outputdir, "2 - D4 dSTORM 1400nm.png"))
	p1 = localizationsplot(r, insetbox = [insetx, insety], include_scalebar = true)
	p2 = insetplot(r, insetx, insety, include_scalebar = true)
	plot(p1, p2, layout = grid(1,2), size=(1024,512), fmt = :png)
end

# ╔═╡ 5db1f140-6ce7-11eb-29f1-ffb93c61c906
let
	# Plot all
	r = localizationsplot.(experimentresults[1][j][i][k] for k ∈ 1:10 for i ∈ 1:4 for j ∈ 1:3)
	p = plot(r..., size=(2048, 2048), layout=grid(10,12), fmt=:png)
	savefig(joinpath(outputdir, "p53-mdm2-all-localizations.png"))
	p
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
	p = [qqnorm(zresid[:, i]) for i ∈ 1:12]
	plot(p..., layout = grid(4, 3), legend = :none, plot_title = "MDM2-p53 qqnorm")
end

# ╔═╡ 6f727ed0-5519-11eb-1462-01c83aa975b6
qqnormplot(mediansflatw)

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
# ╠═f0393540-5514-11eb-2f76-3514cc3bc85c
# ╟─f02df530-5515-11eb-0d9e-8dbb63ed2023
# ╟─1bf98170-5516-11eb-0c24-77ac43672fba
# ╟─25253910-5516-11eb-0a7a-1d73ed51d4a1
# ╟─d8c96c5e-5517-11eb-0861-056d5e425f04
# ╟─87eb8f70-5518-11eb-1639-f91b4bcb5c9b
# ╠═9b7f65a2-5952-11eb-0c21-733df3703bf9
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
# ╠═5f1ce3b0-55c1-11eb-3d03-b53516cfc976
# ╟─756bf2ee-55c1-11eb-2363-d51f3f052b4f
# ╟─b6c33e70-55c1-11eb-3892-d321937a86ed
# ╟─e3220a00-55c1-11eb-080c-451cd462a4d6
# ╟─cbfc9a20-55d0-11eb-3cd3-23a3413a033c
# ╟─2ed8ff90-568e-11eb-10a3-0b97b2425d2f
# ╟─76abe9e0-6af2-11eb-0dc5-4594d89c209b
# ╟─81b932c0-6af2-11eb-0a5c-09c0b3f10519
# ╟─b69b87b2-6af5-11eb-30c2-bfa8159a8e0d
# ╠═deafcdc0-6af4-11eb-1673-0101d277b4b3
# ╟─25d0ce1e-5691-11eb-2cc8-49980b6b2161
# ╟─8c0ebd10-56a9-11eb-3e8c-95b53bf12061
# ╟─faa812c0-56aa-11eb-3b6b-99fe5f200b67
# ╠═012a7c50-56ab-11eb-3d50-c9820925f9ac
# ╟─1637a040-56ac-11eb-2433-49c37383e005
# ╟─87ec3520-56ac-11eb-2a1a-2981155d4168
# ╟─9ce1f5f0-56ac-11eb-22a5-d9b5d6a18978
# ╟─ba7d2760-56ac-11eb-11ca-d543a628c68b
# ╟─c18cc050-56ad-11eb-261c-4177aa31ca7c
# ╠═935aa530-6d73-11eb-2012-4d5e8b594f42
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
# ╠═755d8510-56b3-11eb-3a0c-53b643ba2053
# ╟─bf915350-574e-11eb-2f66-236c4ec9e600
# ╟─d79abcc0-574e-11eb-3777-4f8b85036436
# ╟─f5f23950-574e-11eb-02b5-2fb9b6272ce6
# ╟─27747350-5752-11eb-12b9-dfffc84b15c8
# ╟─248d1560-5753-11eb-2156-89cfb6742205
# ╟─552d3550-6af5-11eb-2fcb-fb755145a298
# ╟─5ab31670-6af5-11eb-2e08-83b5228005d1
# ╟─bc4851c0-6af5-11eb-0db5-cb71d67f93e3
# ╠═af04a0e0-6af5-11eb-05ed-99c27fd1a0ef
# ╟─78883b10-5756-11eb-32cf-b95fdeebda79
# ╟─df7dce20-5756-11eb-3796-037e8e3c49f9
# ╟─faf5aa60-5756-11eb-3b5b-f1ad5e004af2
# ╠═ff845a40-5756-11eb-2ce8-d5ba8458e2ae
# ╟─d2282670-5757-11eb-0bd9-b1aaf0b37d70
# ╟─d99fec80-5757-11eb-39cd-4b0febbf61e3
# ╟─edee22b2-5757-11eb-317a-d3e4d31486cd
# ╟─a05efb90-5758-11eb-10ca-fb9fb8e05821
# ╟─7323e450-5759-11eb-1c62-7ba711bf418b
# ╠═9ecc38b0-6d74-11eb-2d1a-ebfa12f11a64
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
# ╟─b7624d5e-575d-11eb-1221-974537c3b48e
# ╟─8f0dad80-575f-11eb-2312-adae92d50e83
# ╟─a5489880-575f-11eb-1afa-c3bfaf5dccc8
# ╟─4b211122-575f-11eb-0aa9-59fc4e9a2df7
# ╟─408e5d30-5764-11eb-124f-4f2630d1d908
# ╟─d83ab120-6be5-11eb-0b44-176f26137172
# ╠═e5020610-6be5-11eb-25ab-adb9479d6b16
# ╠═840e9e20-6be7-11eb-2b80-37ffe2f4495d
# ╠═7affa0d0-6be8-11eb-3395-555814cab9be
# ╠═13ecfdf0-6bea-11eb-3af4-8d627cccdd2b
# ╠═5db1f140-6ce7-11eb-29f1-ffb93c61c906
# ╟─9c92b290-56ae-11eb-2595-85845fe95f0e
# ╟─2623d480-56af-11eb-2ab8-0d0afd8c45dc
# ╟─a655a120-56ae-11eb-30be-459d7c3067dc
# ╟─7922d1d0-56b0-11eb-261c-5f87ca587394
