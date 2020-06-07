# Simulate a cell:
#    1. Randomly generate positions of two sets of molecules of given amounts m, n
#         - select b = x% of min(m, n) binding pairs, for b pairs, m - b and n - b independent molcules
#         - randomly generate the m - b and n - b molecules
#         - generate the positions of bound pairs d distance apart, then randomly generate a 3D position and squash it
#           to 2D
#             - accounts for antibody distance too
#    2. Randomly generate blinks from each molecule (Poisson?)
#         - accounts for multiple labeling too
#    3. Run the simulated cell through the algorithm and obtain median distance and fraction bound
#    4. Repeat for varying combinations of m, n, x, d
#         - x = { 0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0 }
#         - m = n = { 10, 20, 50, 100, 200 }
#         - d = { 10, 20, 50, 100, 200 } nm

# Basis for choices:
#
# In a particularly dense image, I have observed 2430 molecules and 500 molecules in 25 square microns, and in a less-
# dense area, 301 and 55 molecules, respectively. Note: Using a merge radius of 200 nm practically limits density to
# ~200 molecules in an area this size. Higher density from temporally-displaced localization clouds would reduce
# detection further.
#
# 20 nm  = 10 nm distance + 5 nm nanobodies
# 80 nm  = 70 nm distance + 5 nm nanobodies OR
#          10 nm distance + 35 nm antibody stack
# 140 nm = 70 nm distance + 35 nm antibody stack
#
# For context, microtubulues are 24 nm across, ribosomes 30 nm across
#

fractionsbound = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
moleculecounts = [100, 200, 500, 1000, 2000]
boundradii = [10, 20, 50, 100, 200]

cellradius = 8921 # 250 square micron circle

rootpath = "output/simulation"

using Distributed
using Statistics
using LocalizationMicroscopy
using FileIO
currentworkers = addprocs(exeflags="--project")
@everywhere using SMLMAssociationAnalysis_NCB

function generateunboundmolecules(molecules::Vector{Molecule}, moleculecount, cellradius)
    while length(molecules) < moleculecount
        unboundcoordinates = randomcoordinates2d(moleculecount - length(molecules), cellradius)
        unboundmolecules = [Molecule(Localization(length(molecules) + i, "", unboundcoordinates[1,i], unboundcoordinates[2,i], unboundcoordinates[3,i], 0, 1, 1)) for i ∈ 1:size(unboundcoordinates, 2)]
        append!(molecules, unboundmolecules)
        molecules = merge_close_molecules(molecules, 200)
        foreach(i -> molecules[i].index = i, eachindex(molecules))
    end
    molecules
end

for i ∈ 1:30
    results = ResultSimulate[]

    for moleculecount ∈ moleculecounts
        for fractionbound ∈ fractionsbound
            boundcount = Int(round(fractionbound * moleculecount))

            for boundradius ∈ boundradii
                fractionbound == 0.0 && boundradius != boundradii[1] && break

                println("Processing $moleculecount mols $fractionbound bound $boundradius radius $i")

                boundcoordinates1 = randomcoordinates2d(boundcount, cellradius)
                molecules1 = [Molecule(Localization(i, "", boundcoordinates1[1,i], boundcoordinates1[2,i], boundcoordinates1[3,i], 0, 1, 1)) for i ∈ 1:size(boundcoordinates1, 2)]

                boundcoordinates2 = boundcoordinates1 .+ randomcoordinates2d(boundcount, boundradius)
                molecules2 = [Molecule(Localization(i, "", boundcoordinates2[1,i], boundcoordinates2[2,i], boundcoordinates2[3,i], 0, 1, 1)) for i ∈ 1:size(boundcoordinates2, 2)]

                molecules1 = generateunboundmolecules(molecules1, moleculecount, cellradius)
                molecules2 = generateunboundmolecules(molecules2, moleculecount, cellradius)

                neighbors1, neighbors2, distances = exclusivenearestneighbors(molecules1, molecules2)

                percentileranks = montecarloaffinity(molecules1, molecules2, neighbors1, neighbors2, distances, 200, 10000)
                mediandistance = length(distances) > 0 ? median(distances) : NaN

                data1 = ChannelData("1", molecules1, neighbors1)
                data2 = ChannelData("2", molecules2, neighbors2)
                result = ResultSimulate("", "", 1, "$moleculecount mols $fractionbound bound $boundradius radius", i,
                                [data1, data2], distances, mediandistance, percentileranks)
                push!(results, result)
            end
        end
    end
    save(joinpath(rootpath, "simulationresults$i.jld2"), "results", results)
end

rmprocs(currentworkers)


# plot

nreplicates = 30

using StatsPlots

results = Vector{ResultSimulate}[]
for i ∈ 1:nreplicates
    push!(results, load(joinpath(rootpath, "simulationresults$i.jld2"))["results"])
end

medianmeasurements = Array{Float64,4}(undef, nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))
montecarlomeasurements = Array{Float64,4}(undef, nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))

using InvertedIndices

for i ∈ 1:nreplicates
    replicateresults = results[i]
    lessthanlimit = [(x.distances .< 200) .& (x.percentileranks .< 0.1) for x ∈ replicateresults]
    lessthan10 = count.(lessthanlimit) ./ length.(lessthanlimit)
    mediandistances = map(x -> x.mediandistance, replicateresults)

    medianmeasurements[i,:,1,:] .= reshape(repeat(mediandistances[1:36:end], inner=length(boundradii)), length(boundradii), length(moleculecounts))
    medianmeasurements[i,:,2:end,:] .= reshape(mediandistances[Not(1:36:end)], length(boundradii), length(fractionsbound)-1, length(moleculecounts))

    montecarlomeasurements[i,:,1,:] .= reshape(repeat(lessthan10[1:36:end], inner=length(boundradii)), length(boundradii), length(moleculecounts))
    montecarlomeasurements[i,:,2:end,:] .= reshape(lessthan10[Not(1:36:end)], length(boundradii), length(fractionsbound)-1, length(moleculecounts))
end

import Plots.mm

# monte carlo plot
p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
xpos = reshape(repeat(["0", "1", "2", "5", "10", "20", "50", "100"], inner=5 * nreplicates, outer=5), nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        p[i,j] = boxplot(xpos[:,j,:,i], montecarlomeasurements[:,j,:,i], legend=:none, outliers=false, tickfontsize=12, seriescolor=:white, yaxis=(0:1), xrotation=60, left_margin=-0.5mm, bottom_margin=-3mm, gridopacity=0.2)
        i > 1 && plot!(p[i,j], yaxis = false, left_margin=-10mm)
        j < length(boundradii) && plot!(p[i,j], xaxis = false, bottom_margin=-10mm)
        dotplot!(p[i,j], xpos[:,j,:,i], montecarlomeasurements[:,j,:,i], legend=:none, bar_width=0.1, markerstrokewidth=0, markersize=2, markercolor=:black)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(1024,1024))
savefig(joinpath(rootpath, "simulation_montecarlo.png"))

# median plot
p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        p[i,j] = boxplot(xpos[:,j,:,i], medianmeasurements[:,j,:,i], legend=:none, outliers=false, yaxis=((0, 1200), 0:200:1200), tickfontsize=12, seriescolor=:white, xrotation=60, left_margin=-0.5mm, bottom_margin=-3mm, gridopacity=0.2)
        i > 1 && plot!(p[i,j], yaxis = false, left_margin=-10mm)
        j < length(boundradii) && plot!(p[i,j], xaxis = false, bottom_margin=-10mm)
        dotplot!(p[i,j], xpos[:,j,:,i], medianmeasurements[:,j,:,i], legend=:none, bar_width=0.1, markerstrokewidth=0, markersize=2, markercolor=:black)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(1024,1024))
savefig(joinpath(rootpath,"simulation_median.png"))

# normalized monte carlo plot

zero_montecarlomeasurements = reshape(repeat(median(montecarlomeasurements[:,:,1,:], dims = 1), inner = (30, 8, 1)), 30, 5, 8, 5)
one_montecarlomeasurements = reshape(repeat(median(montecarlomeasurements[:,:,length(fractionsbound),:], dims = 1), inner = (30, 1, 1), outer = (1, 8, 1)), 30, 5, 8, 5)
normalizedmontecarlomeasurements = (montecarlomeasurements .- zero_montecarlomeasurements) ./ (one_montecarlomeasurements .- zero_montecarlomeasurements)

p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
xpos = reshape(repeat(["0", "1", "2", "5", "10", "20", "50", "100"], inner=5 * nreplicates, outer=5), nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        p[i,j] = boxplot(xpos[:,j,:,i], normalizedmontecarlomeasurements[:,j,:,i], legend=:none, outliers=false, tickfontsize=12, seriescolor=:white, yaxis=(0:1), xrotation=60, left_margin=-0.5mm, bottom_margin=-3mm, gridopacity=0.2)
        plot!(p[i,j], xpos[1,j,:,i], fractionsbound, lineopacity = 0.25, linecolor=:black, linewidth=3)
        i > 1 && plot!(p[i,j], yaxis = false, left_margin=-10mm)
        j < length(boundradii) && plot!(p[i,j], xaxis = false, bottom_margin=-10mm)
        dotplot!(p[i,j], xpos[:,j,:,i], normalizedmontecarlomeasurements[:,j,:,i], legend=:none, bar_width=0.1, markerstrokewidth=0, markersize=2, markercolor=:black)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(1024,1024))
savefig(joinpath(rootpath,"simulation_montecarlo_normalized.png"))



# normalized monte carlo plot deviation

deviations = (median(normalizedmontecarlomeasurements[:,:,eachindex(fractionsbound)[Not(1)],:], dims = 1) .- reshape(repeat(fractionsbound[Not(1)], inner = 5, outer = 5), 1, 5, 7, 5)) ./ reshape(repeat(fractionsbound[Not(1)], inner = 5, outer = 5), 1, 5, 7, 5)

p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
xpos = reshape(repeat(["1", "2", "5", "10", "20", "50", "100"], inner=5, outer=5), length(boundradii), length(fractionsbound) - 1, length(moleculecounts))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        p[i,j] = plot(xpos[j,:,i], deviations[1,j,:,i], legend=:none, outliers=false, tickfontsize=12, yaxis=((-5, 2), -5:1:2), xrotation=60, left_margin=-0.5mm, bottom_margin=-3mm, gridopacity=0.2, linewidth=3, linecolor=:black)
        i > 1 && plot!(p[i,j], yaxis = false, left_margin=-10mm)
        j < length(boundradii) && plot!(p[i,j], xaxis = false, bottom_margin=-10mm)
        plot!(p[i,j], xpos[j,:,i], repeat([0], 7), linecolor=:gray, linewidth=3, lineopacity=0.5)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(1024,1024))
savefig(joinpath(rootpath,"simulation_montecarlo_normalized_deviations.png"))
