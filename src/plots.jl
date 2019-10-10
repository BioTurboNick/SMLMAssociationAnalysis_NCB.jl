import Plots.scatter
import Plots.mm

function showPairs(molecules1, molecules2, pairedmolecules1, pairedmolecules2, distances)
    if length(pairedmolecules1) != length(pairedmolecules2) != length(distances)
        throw(ErrorException("Lengths of paired inputs must be identical."))
    end

    p = plot()
    pairedcoordinates1 = Main.extractcoordinates(pairedmolecules1)
    pairedcoordinates2 = Main.extractcoordinates(pairedmolecules2)
    coordinates1 = Main.extractcoordinates(molecules1)
    coordinates2 = Main.extractcoordinates(molecules2)
    scatter!(coordinates1[1,:], coordinates1[2,:], markersize = 2, markercolor = :white)
    scatter!(coordinates2[1,:], coordinates2[2,:], markersize = 2, markercolor = :white)
    scatter!(pairedcoordinates1[1,:], pairedcoordinates1[2,:], markersize = 1, markerstrokecolor=nothing, markercolor = :red)
    scatter!(pairedcoordinates2[1,:], pairedcoordinates2[2,:], markersize = 1, markerstrokecolor=nothing, markercolor = :green)
end

using Printf

function distanceprobabilityplot(result::Result)
    lessthan10 = count((result.percentileranks .< 0.1) .& (result.distances .< 200)) / length(result.percentileranks)
    p1 = scatter(result.distances, result.percentileranks,
                 xaxis=("distance (nm)", (0,1000), 45),
                 yaxis=(:hide),
                 marker=(2, stroke(0)),
                 legend=:none,
                 left_margin=-30mm,
                 top_margin=-30mm,
                 annotations=(190,0,text(Printf.@sprintf("%.3f", lessthan10), 8, :right, :bottom)),
                 ywiden=false,
                 tick_direction=:out)
    line = plot!(p1, [0; 1000], [0.1, 0.1], line=(:gray))
    line = plot!(p1, [200; 200], [0, 1.0], line=(:gray))
    p2 = scatter(result.distances, result.percentileranks,
                 xaxis=((1000,15000), (1000:7000:15000), 45),
                 xformatter=:plain,
                 yaxis=(:hide),
                 marker=(2, stroke(0)),
                 legend=:none,
                 top_margin=-30mm,
                 left_margin=-20mm,
                 ywiden=false,
                 tick_direction=:out)

    p3 = histogram(result.distances,
                   bins=0:25:1000,
                   left_margin=-10mm,
                   xaxis=((0,1000), :hide),
                   linealpha=0,
                   legend=:none,
                   tick_direction=:out)
    p4 = histogram(result.distances,
                   bins=1000:25:15000,
                   xaxis=((1000,15000), (1000:7000:15000), :hide),
                   legend=:none,
                   linealpha=0,
                   left_margin=-30mm,
                   yaxis=(:hide))
    p5 = histogram(result.percentileranks,
                   bins=0:0.05:1.01,
                   orientation = :horizontal,
                   linealpha=0,
                   xaxis=(:flip, 45),
                   yaxis=("p(chance association)"),
                   legend=:none,
                   top_margin=-30mm,
                   ywiden=false,
                   tick_direction=:out)
    plot(Plots.Plot(), p3, p4, p5, p1, p2, layout = grid(2, 3, widths=[0.15,0.65,0.15], heights=[0.1,0.9]), link=:y)
end

using LocalizationMicroscopy
function moleculesplot(result::Result)
    mol1coords = result.channels[1].molecules |> LocalizationMicroscopy.extractcoordinates
    mol2coords = result.channels[2].molecules |> LocalizationMicroscopy.extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(2, stroke(0), :red))
    plot!([1000; 6000], [1000; 1000], line=(3, :black), annotations=(1000,1250,text("5 \\mum", 10, :left, :top)))

    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(2, stroke(0), :green))
    plot!(aspect_ratio=:equal, xlims=(0, 40960), yaxis=((0, 40960), :flip), legend=:none, grid=:hide, ticks=(0),
          framestyle=:box) # check limit
end

function localizationsplot(result::Result; insetbox = [[0,0], [0,0]])
    loc1coords = mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules) |> extractcoordinates
    loc2coords = mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules) |> extractcoordinates
    scatter(loc1coords[1,:], loc1coords[2,:], marker=(1, stroke(0), 0.75, :red), framestyle=:none)
    plot!([1000; 6000], [1000; 1000], line=(3, :black), annotations=(1000,1250,text("5 \\mum", 10, :left, :top)))

    scatter!(loc2coords[1,:], loc2coords[2,:], marker=(1, stroke(0), 0.75, :green))
    plot!(aspect_ratio=:equal, xlims=(0, 40960), yaxis=((0, 40960), :flip), legend=:none, grid=:hide, ticks=(0))

    if insetbox != [[0,0], [0,0]]
        xlims, ylims = first(insetbox), last(insetbox)
        plot!(xlims, repeat([first(ylims)], 2), line=(1, :black))
        plot!(xlims, repeat([last(ylims)], 2), line=(1, :black))
        plot!(repeat([first(xlims)], 2), ylims, line=(1, :black))
        plot!(repeat([last(xlims)], 2), ylims, line=(1, :black))
    end

    plot!()
end

function localizationsplot_forprint(result::Result; color1 = :red, color2 = :green, insetbox = [[0,0], [0,0]])
    loc1coords = mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules) |> extractcoordinates
    loc2coords = mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules) |> extractcoordinates
    scatter(loc1coords[1,:], loc1coords[2,:], marker=(8, stroke(0), 0.75, color1), size=(2048,2048), framestyle=:none)
    plot!([1000; 6000], [1000; 1000], line=(3, :black), annotations=(1000,1250,text("5 \\mum", 10, :left, :top)))

    scatter!(loc2coords[1,:], loc2coords[2,:], marker=(8, stroke(0), 0.75, color2))
    plot!(aspect_ratio=:equal, xlims=(0, 40960), yaxis=((0, 40960), :flip), legend=:none, grid=:hide, ticks=(0))

    if insetbox != [[0,0], [0,0]]
        xlims, ylims = first(insetbox), last(insetbox)
        plot!(xlims, repeat([first(ylims)], 2), line=(8, :black))
        plot!(xlims, repeat([last(ylims)], 2), line=(8, :black))
        plot!(repeat([first(xlims)], 2), ylims, line=(8, :black))
        plot!(repeat([last(xlims)], 2), ylims, line=(8, :black))
    end

    plot!()
end

ENV["GKS_ENCODING"] = "utf-8"

function moleculesinsetplot(result::Result, xlims, ylims)
    mol1coords = result.channels[1].molecules |> LocalizationMicroscopy.extractcoordinates
    mol2coords = result.channels[2].molecules |> LocalizationMicroscopy.extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(4, stroke(0), :red))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(4, stroke(0), :green))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0), framestyle=:box)
end

function moleculesinsetplot_forprint(result::Result, xlims, ylims; color1 = :red, color2 = :green)
    mol1coords = result.channels[1].molecules |> LocalizationMicroscopy.extractcoordinates
    mol2coords = result.channels[2].molecules |> LocalizationMicroscopy.extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(16, stroke(0), color1), size=(2048,2048), framestyle=:none)
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(16, stroke(0), color2))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0))
end

function localizationsinsetplot(result::Result, xlims, ylims)
    loc1coords = mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules) |> extractcoordinates
    loc2coords = mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules) |> extractcoordinates
    scatter(loc1coords[1,:], loc1coords[2,:], marker=(1, stroke(0), :red))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(loc2coords[1,:], loc2coords[2,:], marker=(1, stroke(0), :green))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0), framestyle=:box)
end

function localizationsinsetplot_forprint(result::Result, xlims, ylims; color1 = :red, color2 = :green)
    loc1coords = mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules) |> extractcoordinates
    loc2coords = mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules) |> extractcoordinates
    scatter(loc1coords[1,:], loc1coords[2,:], marker=(8, stroke(0), 0.75, color1), framestyle=:none, size=(2048,2048))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(loc2coords[1,:], loc2coords[2,:], marker=(8, stroke(0), 0.75, color2))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0))
end

function insetplot(result::Result, xlims, ylims)
    mol1coords = result.channels[1].molecules |> LocalizationMicroscopy.extractcoordinates
    mol2coords = result.channels[2].molecules |> LocalizationMicroscopy.extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(8, stroke(2, :red), :white))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(8, stroke(2, :green), :white))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0))

    loc1coords = mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules) |> extractcoordinates
    loc2coords = mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules) |> extractcoordinates
    scatter!(loc1coords[1,:], loc1coords[2,:], marker=(1, stroke(0), :red))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(loc2coords[1,:], loc2coords[2,:], marker=(1, stroke(0), :green))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0), framestyle=:box)
end

function insetplot_forprint(result::Result, xlims, ylims; color1 = :red, color2 = :green)
    mol1coords = result.channels[1].molecules |> LocalizationMicroscopy.extractcoordinates
    mol2coords = result.channels[2].molecules |> LocalizationMicroscopy.extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(24, stroke(2, color1), 0.75, :white),
            framestyle=:none, size=(2048,2048))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(24, stroke(2, color2), 0.75, :white))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0))

    loc1coords = mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules) |> extractcoordinates
    loc2coords = mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules) |> extractcoordinates
    scatter!(loc1coords[1,:], loc1coords[2,:], marker=(8, stroke(0), 0.75, color1))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(loc2coords[1,:], loc2coords[2,:], marker=(8, stroke(0), 0.75, color2))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0))
end

function yzframelocalizationplot(result::Result)
    loc1coords = mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules) |> extractcoordinates
    loc1frames = [y.frame for y ∈ mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules)]
    loc2coords = mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules) |> extractcoordinates
    loc2frames = [y.frame for y ∈ mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules)]

    scatter(loc1frames, loc1coords[2,:], marker=(2, :red, stroke(0)), xaxis=("frames", (0,22000)),
            yaxis=("y", (0,40960), :flip, 90), legend=:none)
    scatter!(loc2frames, loc2coords[2,:], marker=(2, :green, stroke(0)))
    plot!(tick_direction=:out, size=(1024,512))
end
