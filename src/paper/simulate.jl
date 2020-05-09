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

#cellradius = 2821 # 25 square micron circle
#cellradius = 3989 # 50 square micron circle
#cellradius = 5642 # 100 square micron circle
cellradius = 8921 # 250 square micron circle

rootpath = "C:/Users/nicho/Dropbox (Partners HealthCare)/STORM MATLAB/STORM Single Molecule Clustering/MonteCarloAffinity/Simulated data"

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
    results = Result[]

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

                percentileranks = montecarloaffinity(molecules1, molecules2, neighbors1, neighbors2, distances, 800, 10000)
                mediandistance = length(distances) > 0 ? median(distances) : NaN

                data1 = ChannelData("1", molecules1, neighbors1)
                data2 = ChannelData("2", molecules2, neighbors2)
                result = Result("", "", 1, "$moleculecount mols $fractionbound bound $boundradius radius", i,
                                [data1, data2], distances, mediandistance, percentileranks)
                push!(results, result)
            end
        end
    end
    save(joinpath(rootpath, "simulationresults_new$i.jld2"), "results", results)
end

rmprocs(currentworkers)


# plot

nreplicates = 30

using StatsPlots

results = Vector{Result}[]
for i ∈ 1:nreplicates
    push!(results, load(joinpath(rootpath, "simulationresults_new$i.jld2"))["results"])
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
p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
xpos = reshape(repeat(["0", "1", "2", "5", "10", "20", "50", "100"], inner=5 * nreplicates, outer=5), nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        p[i,j] = boxplot(xpos[:,j,:,i], montecarlomeasurements[:,j,:,i], legend=:none, outliers=false, seriescolor=:white, linecolor=:gray, yaxis=((0.0,1.0), 0:0.2:1), tickfontsize=18, left_margin=5mm, xrotation=45)
        i > 1 && plot!(ytickfontcolor=:white)
        j < 5 && plot!(xtickfontcolor=:white)
        dotplot!(p[i,j], xpos[:,j,:,i], montecarlomeasurements[:,j,:,i], legend=:none, markerstrokewidth=0, markersize=3, markercolor=:black)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(2048,2048))
savefig("simulation_montecarlo.png")

p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        p[i,j] = boxplot(xpos[:,j,:,i], medianmeasurements[:,j,:,i], legend=:none, outliers=false, seriescolor=:white, linecolor=:gray, yaxis=((0, 1200), 0:200:1200), tickfontsize=18, left_margin=10mm, xrotation=45)
        i > 1 && plot!(ytickfontcolor=:white)
        j < 5 && plot!(xtickfontcolor=:white)
        dotplot!(p[i,j], xpos[:,j,:,i], medianmeasurements[:,j,:,i], legend=:none, markerstrokewidth=0, markersize=3, markercolor=:black)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(2048,2048))
savefig("simulation_median.png")


fractionsboundmatrix = reshape(repeat(fractionsbound, inner=5 * nreplicates, outer=5), nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))
montecarlodeviations = fractionsboundmatrix .- montecarlomeasurements

p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        p[i,j] = boxplot(xpos[:,j,:,i], montecarlodeviations[:,j,:,i], legend=:none, outliers=false, seriescolor=:white, linecolor=:gray, yaxis=((-0.1,1), -0.1:0.1:1), tickfontsize=18, left_margin=10mm, xrotation=45)
        i > 1 && plot!(ytickfontcolor=:white)
        j < 5 && plot!(xtickfontcolor=:white)
        dotplot!(p[i,j], xpos[:,j,:,i], montecarlodeviations[:,j,:,i], legend=:none, markerstrokewidth=0, markersize=3, markercolor=:black)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(2048,2048))
savefig("simulation_montecarlodeviation.png")

montecarlorelativedeviations = montecarlodeviations[:,:,2:8,:] ./ fractionsboundmatrix[:,:,2:8,:]

p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        p[i,j] = boxplot(xpos[:,j,2:8,i], montecarlorelativedeviations[:,j,:,i], legend=:none, outliers=false, seriescolor=:white, linecolor=:gray, yaxis=((-10,2), -10:2:2), tickfontsize=18, left_margin=6mm, xrotation=45)
        i > 1 && plot!(ytickfontcolor=:white)
        j < 5 && plot!(xtickfontcolor=:white)
        dotplot!(p[i,j], xpos[:,j,2:8,i], montecarlorelativedeviations[:,j,:,i], legend=:none, markerstrokewidth=0, markersize=3, markercolor=:black)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(2048,2048))
savefig("simulation_montecarlorelativedeviation.png")

using Distributions
dists = Normal.(mean(montecarlomeasurements, dims=1), std(montecarlomeasurements, dims=1))

boundradiimatrix = reshape(repeat(boundradii, inner=nreplicates, outer=5*8), nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))
moleculedensitymatrix = reshape(repeat(moleculecounts ./ 250, inner=nreplicates * 5, outer=8), nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))

using DataFrames
df = DataFrame(X = fractionsboundmatrix |> vec, Y = montecarlomeasurements |> vec, R = boundradiimatrix |> vec, D = moleculedensitymatrix |> vec)
ols = lm(@formula(Y ~ (X^2) / (D * R)), df)
predicted = DataFrame(X = fractionsboundmatrix[1,:,:,:] |> vec, R = boundradiimatrix[1,:,:,:] |> vec, D = moleculedensitymatrix[1,:,:,:] |> vec)
predicted.Y = predict(ols, predicted)

p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        p[i,j] = boxplot(xpos[:,j,:,i], montecarlomeasurements[:,j,:,i], legend=:none, outliers=false, seriescolor=:white, linecolor=:gray, yaxis=((0.0,1.0), 0:0.2:1), tickfontsize=18, left_margin=5mm, xrotation=45)
        i > 1 && plot!(ytickfontcolor=:white)
        j < 5 && plot!(xtickfontcolor=:white)
        dotplot!(p[i,j], xpos[:,j,:,i], montecarlomeasurements[:,j,:,i], legend=:none, markerstrokewidth=0, markersize=3, markercolor=:black)
        plot!(p[i,j], xpos[1,j,:,i], reshape(predicted.Y, length(boundradii), length(fractionsbound), length(moleculecounts))[j,:,i])
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(2048,2048))
savefig("simulation_montecarlo_regression.png")


p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
xpos = reshape(repeat([0, .01, .02, .05, .1, .2, .5, 1], inner=5 * nreplicates, outer=5), nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        df = DataFrame(X = fractionsboundmatrix[:,j,:,i] |> vec, Y = montecarlomeasurements[:,j,:,i] |> vec, R = boundradiimatrix[:,j,:,i] |> vec, D = moleculedensitymatrix[:,j,:,i] |> vec)
        ols = lm(@formula(Y ~ X), df)
        predicted = DataFrame(X = fractionsboundmatrix[1,j,:,i] |> vec, D = moleculedensitymatrix[1,j,:,i] |> vec) #, R = boundradiimatrix[1,j,:,i] |> vec, )
        predicted.Y = predict(ols, predicted)
        p[i,j] = boxplot(xpos[:,j,:,i], montecarlomeasurements[:,j,:,i], legend=:none, outliers=false, seriescolor=:white, bar_width=0.1, linecolor=:gray, yaxis=((0.0,1.0), 0:0.2:1), tickfontsize=18, left_margin=5mm, xrotation=45)
        i > 1 && plot!(ytickfontcolor=:white)
        j < 5 && plot!(xtickfontcolor=:white)
        dotplot!(p[i,j], xpos[:,j,:,i], montecarlomeasurements[:,j,:,i], legend=:none, markerstrokewidth=0, markersize=3, bar_width=0.1, markercolor=:black)
        plot!(p[i,j], xpos[1,j,:,i], predicted.Y)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(2048,2048))
savefig("simulation_montecarlo_regression.png")


p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
xpos = reshape(repeat(["0", "1", "2", "5", "10", "20"], inner=5 * nreplicates, outer=5), nreplicates, length(boundradii), length(fractionsbound) - 2, length(moleculecounts))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        df = DataFrame(X = fractionsboundmatrix[:,j,:,i] |> vec, Y = montecarlomeasurements[:,j,:,i] |> vec, R = boundradiimatrix[:,j,:,i] |> vec, D = moleculedensitymatrix[:,j,:,i] |> vec)
        ols = lm(@formula(Y ~ exp(X) + X), df)
        predicted = DataFrame(X = fractionsboundmatrix[1,j,1:6,i] |> vec, D = moleculedensitymatrix[1,j,1:6,i] |> vec) #, R = boundradiimatrix[1,j,:,i] |> vec, )
        predicted.Y = predict(ols, predicted)
        p[i,j] = boxplot(xpos[:,j,:,i], montecarlomeasurements[:,j,1:6,i], legend=:none, outliers=false, seriescolor=:white, linecolor=:gray, yaxis=((0.0,0.3), 0:0.05:0.3), tickfontsize=18, left_margin=10mm, xrotation=45)
        i > 1 && plot!(ytickfontcolor=:white)
        j < 5 && plot!(xtickfontcolor=:white)
        dotplot!(p[i,j], xpos[:,j,:,i], montecarlomeasurements[:,j,1:6,i], legend=:none, markerstrokewidth=0, markersize=3, markercolor=:black)
        plot!(p[i,j], xpos[1,j,:,i], predicted.Y)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(2048,2048))
savefig("simulation_montecarlo_0_20_regression.png")

p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
xpos = reshape(repeat(["0", "1", "2", "5", "10", "20", "50", "100"], inner=5 * nreplicates, outer=5), nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        minval = minimum(mean(montecarlomeasurements[:,j,:,i], dims = 1))
        normalized = (montecarlomeasurements[:,j,:,i] .- minval) ./ (mean(montecarlomeasurements[:,j,8,i]) - minval)
        p[i,j] = boxplot(xpos[:,j,:,i], normalized, legend=:none, outliers=false, seriescolor=:white, linecolor=:gray, yaxis=((0.0,1.0), 0:0.2:1), tickfontsize=18, left_margin=5mm, xrotation=45)
        i > 1 && plot!(ytickfontcolor=:white)
        j < 5 && plot!(xtickfontcolor=:white)
        dotplot!(p[i,j], xpos[:,j,:,i], normalized, legend=:none, markerstrokewidth=0, markersize=3, markercolor=:black)
        plot!(p[i,j], xpos[1,j,:,i], fractionsboundmatrix[1,j,1:8,i], opacity=0.5, linewidth=10)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(2048,2048))
savefig("simulation_montecarlo_normalized.png")

p = Array{Plots.Plot,2}(undef, length(moleculecounts), length(boundradii))
xpos = reshape(repeat(["0", "1", "2", "5", "10", "20", "50", "100"], inner=5 * nreplicates, outer=5), nreplicates, length(boundradii), length(fractionsbound), length(moleculecounts))
for i ∈ eachindex(moleculecounts)
    for j ∈ eachindex(boundradii)
        minval = minimum(mean(montecarlomeasurements[:,j,:,i], dims = 1))
        normalized = median((montecarlomeasurements[:,j,:,i] .- minval) ./ (mean(montecarlomeasurements[:,j,8,i]) - minval), dims = 1)
        p[i,j] = plot(xpos[1,1,2:8,1], repeat([0], 7), linecolor=:black, linewidth=1)
        plot!(xpos[1,j,2:8,i], (normalized[2:8] .- fractionsboundmatrix[1,j,2:8,i]) ./ fractionsboundmatrix[1,j,2:8,i], legend=:none, linewidth=10, seriescolor=:white, linecolor=:gray, yaxis=((-2, 2), -2:0.5:2), tickfontsize=18, left_margin=5mm, xrotation=45)
        i > 1 && plot!(ytickfontcolor=:white)
        j < 5 && plot!(xtickfontcolor=:white)
    end
end

plot(p..., layout=grid(length(moleculecounts), length(boundradii)), size=(2048,2048))
savefig("simulation_montecarlo_normalized_deviations.png")
