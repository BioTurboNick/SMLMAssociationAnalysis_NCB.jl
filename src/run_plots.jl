using SMLMAssociationAnalysis_NCB

# PARAMETERS

outputpath = "path/to/your/output/results.jld2"

percentilerankthreshold = 0.1
distancethreshold = 200

# END PARAMETERS

results = FileIO.load(outputpath)["results"]

# generate localization plots (not print quality)
for result ∈ results
    p = localizationsplot(results)
    savefig(p, "$(result.cell) localizations.png")
end

# generate molecule plots (not print quality)
for result ∈ results
    p = moleculesplot(results)
    savefig(p, "$(result.cell) molecules.png")
end

# generate distance-probability plots
for result ∈ results
    p = distanceprobabilityplot(results)
    savefig(p, "$(result.cell) distance-probability.png")
end

# calculate cell summary data (median distance, fraction bound)
percentilerankslessthan10replicate = [x.percentileranks .< percentilerankthreshold for x ∈ results]
bound = [(x.distances .< distancethreshold) .& (x.percentileranks .< percentilerankthreshold) for x ∈ replicateresults]
fractionbound = count.(bound) ./ length.(percentilerankslessthan10replicate)
mediandistances = map(x -> x.mediandistance, replicateresults)

localization1counts = [map(x -> x.group.localizations |> length, y.channels[1].molecules) |> sum for y ∈ replicateresults]
localization2counts = [map(x -> x.group.localizations |> length, y.channels[2].molecules) |> sum for y ∈ replicateresults]
molecule1counts = map(x -> x.channels[1].molecules |> length, replicateresults)
molecule2counts = map(x -> x.channels[2].molecules |> length, replicateresults)
