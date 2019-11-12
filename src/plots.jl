import Plots.scatter
import Plots.mm

using LocalizationMicroscopy

function showpairs(molecules1, molecules2, pairedmolecules1, pairedmolecules2, distances)
    if length(pairedmolecules1) != length(pairedmolecules2) != length(distances)
        throw(ErrorException("Lengths of paired inputs must be identical."))
    end

    p = plot()
    pairedcoordinates1 = extractcoordinates(pairedmolecules1)
    pairedcoordinates2 = extractcoordinates(pairedmolecules2)
    coordinates1 = extractcoordinates(molecules1)
    coordinates2 = extractcoordinates(molecules2)
    scatter!(coordinates1[1,:], coordinates1[2,:], markersize = 2, markercolor = :white)
    scatter!(coordinates2[1,:], coordinates2[2,:], markersize = 2, markercolor = :white)
    scatter!(pairedcoordinates1[1,:], pairedcoordinates1[2,:], markersize = 1, markerstrokecolor=nothing, markercolor = :magenta)
    scatter!(pairedcoordinates2[1,:], pairedcoordinates2[2,:], markersize = 1, markerstrokecolor=nothing, markercolor = :green)
end

using Printf

distanceprobabilityplot(result::Result) = distanceprobabilityplot(result.distances, result.percentileranks)
function distanceprobabilityplot(distances::Vector{Float64}, percentileranks::Vector{Float64})
    lessthan10 = count((percentileranks .< 0.1) .& (distances .< 200)) / length(percentileranks)
    p1 = scatter(distances, percentileranks,
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
    p2 = scatter(distances, percentileranks,
                 xaxis=((1000,15000), (1000:7000:15000), 45),
                 xformatter=:plain,
                 yaxis=(:hide),
                 marker=(2, stroke(0)),
                 legend=:none,
                 top_margin=-30mm,
                 left_margin=-20mm,
                 ywiden=false,
                 tick_direction=:out)

    p3 = histogram(distances,
                   bins=0:25:1000,
                   left_margin=-10mm,
                   xaxis=((0,1000), :hide),
                   linealpha=0,
                   legend=:none,
                   tick_direction=:out)
    p4 = histogram(distances,
                   bins=1000:25:15000,
                   xaxis=((1000,15000), (1000:7000:15000), :hide),
                   legend=:none,
                   linealpha=0,
                   left_margin=-30mm,
                   yaxis=(:hide))
    p5 = histogram(percentileranks,
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

moleculesplot(result::Result) = moleculesplot(result.channels[1].molecules, result.channels[2].molecules)
function moleculesplot(molecules1::Vector{Molecule}, molecules2::Vector{Molecule})
    mol1coords = molecules1 |> extractcoordinates
    mol2coords = molecules2 |> extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(2, stroke(0), :magenta))
    plot!([1000; 6000], [1000; 1000], line=(3, :black), annotations=(1000,1250,text("5 \\mum", 10, :left, :top)))

    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(2, stroke(0), :green))
    plot!(aspect_ratio=:equal, xlims=(0, 40960), yaxis=((0, 40960), :flip), legend=:none, grid=:hide, ticks=(0),
          framestyle=:box) # check limit
end

neighborsplot(result::Result) = neighborsplot(result.channels[1].molecules, result.channels[2].molecules, result.channels[1].neighbormolecules, result.channels[2].neighbormolecules, result.percentileranks, result.distances)
function neighborsplot(molecules1::Vector{Molecule}, molecules2::Vector{Molecule}, neighbormolecules1::Vector{Molecule}, neighbormolecules2::Vector{Molecule}, percentileranks::Vector{AbstractFloat}, distances::Vector{AbstractFloat})
    if length(neighbormolecules1) != length(neighbormolecules2) != length(percentileranks) != length(distances)
        throw(ErrorException("Lengths of paired inputs must be identical."))
    end
    mol1coords = molecules1 |> extractcoordinates
    mol2coords = molecules2 |> extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(2, stroke(0), :magenta, 0.25))
    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(2, stroke(0), :green, 0.25))

    mol1coords = neighbormolecules1 |> extractcoordinates
    mol2coords = neighbormolecules2 |> extractcoordinates
    scatter!(mol1coords[1,:], mol1coords[2,:], marker=(2, stroke(0), :magenta, 0.25))
    plot!([1000; 6000], [1000; 1000], line=(3, :black), annotations=(1000,1250,text("5 \\mum", 10, :left, :top)))
    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(2, stroke(0), :green, 0.25))

    isbound = (percentileranks .< 0.1) .& (distances .< 200)
    boundmolecules1 = neighbormolecules1[isbound]
    boundmolecules2 = neighbormolecules2[isbound]
    mol1coords = boundmolecules1 |> extractcoordinates
    mol2coords = boundmolecules2 |> extractcoordinates
    scatter!(mol1coords[1,:], mol1coords[2,:], marker=(2, stroke(0), :black, 0.5))
    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(2, stroke(0), :black, 0.5))

    plot!(aspect_ratio=:equal, xlims=(0, 40960), yaxis=((0, 40960), :flip), legend=:none, grid=:hide, ticks=(0),
          framestyle=:box) # check limit
end

neighborsplot_forprint(result::Result, include_scalebar=true) = neighborsplot_forprint(result.channels[1].molecules, result.channels[2].molecules, result.channels[1].neighbormolecules, result.channels[2].neighbormolecules, result.percentileranks, result.distances, true)
function neighborsplot_forprint(molecules1::Vector{Molecule}, molecules2::Vector{Molecule}, neighbormolecules1::Vector{Molecule}, neighbormolecules2::Vector{Molecule}, percentileranks::Vector{AbstractFloat}, distances::Vector{AbstractFloat}, include_scalebar=true)
    if length(neighbormolecules1) != length(neighbormolecules2) != length(percentileranks) != length(distances)
        throw(ErrorException("Lengths of paired inputs must be identical."))
    end
    mol1coords = molecules1 |> extractcoordinates
    mol2coords = molecules2 |> extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(8, stroke(0), :magenta, 0.25))
    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(8, stroke(0), :green, 0.25))
    mol1coords = neighbormolecules1 |> extractcoordinates
    mol2coords = neighbormolecules2 |> extractcoordinates

    scatter!(mol1coords[1,:], mol1coords[2,:], marker=(8, stroke(0), :magenta, 0.25))
    plot!([1000; 6000], [1000; 1000], line=(3, :black), annotations=(1000,1250,text("5 \\mum", 10, :left, :top)))
    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(8, stroke(0), :green, 0.25))

    isbound = (percentileranks .< 0.1) .& (distances .< 200)
    boundmolecules1 = neighbormolecules1[isbound]
    boundmolecules2 = neighbormolecules2[isbound]
    mol1coords = boundmolecules1 |> extractcoordinates
    mol2coords = boundmolecules2 |> extractcoordinates
    scatter!(mol1coords[1,:], mol1coords[2,:], marker=(8, stroke(0), :black, 0.5))
    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(8, stroke(0), :black, 0.5))

    plot!(aspect_ratio=:equal, xlims=(0, 40960), yaxis=((0, 40960), :flip), legend=:none, grid=:hide, ticks=(0),
          framestyle=:box, size=(2048,2048)) # check limit
end

localizationsplot(result::Result; kwargs...) = localizationsplot(result.channels[1].molecules, result.channels[2].molecules; kwargs...)
localizationsplot(molecules1, molecules2; kwargs...) = localizationsplot(mapreduce(x -> x.group.localizations, vcat, molecules1), mapreduce(x -> x.group.localizations, vcat, molecules2); kwargs...)
function localizationsplot(localizations1::Vector{Localization}, localizations2::Vector{Localization}; insetbox = [[0,0], [0,0]])
    loc1coords = localizations1 |> extractcoordinates
    loc2coords = localizations2 |> extractcoordinates
    scatter(loc1coords[1,:], loc1coords[2,:], marker=(1, stroke(0), 0.75, :magenta), framestyle=:none)
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

function localizationsplot_forprint(result::Result; color1 = :magenta, color2 = :green, insetbox = [[0,0], [0,0]], include_scalebar = false)
    loc1coords = mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules) |> extractcoordinates
    loc2coords = mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules) |> extractcoordinates
    scatter(loc1coords[1,:], loc1coords[2,:], marker=(8, stroke(0), 0.75, color1), size=(2048,2048), framestyle=:none)
    include_scalebar && plot!([1000; 6000], [1000; 1000], line=(3, :black), annotations=(1000,1250,text("5 \\mum", 10, :left, :top)))

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
    mol1coords = result.channels[1].molecules |> extractcoordinates
    mol2coords = result.channels[2].molecules |> extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(4, stroke(0), :magenta))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(4, stroke(0), :green))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0), framestyle=:box)
end

function moleculesinsetplot_forprint(result::Result, xlims, ylims; color1 = :magenta, color2 = :green)
    mol1coords = result.channels[1].molecules |> extractcoordinates
    mol2coords = result.channels[2].molecules |> extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(16, stroke(0), color1), size=(2048,2048), framestyle=:none)
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(16, stroke(0), color2))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0))
end

function localizationsinsetplot(result::Result, xlims, ylims)
    loc1coords = mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules) |> extractcoordinates
    loc2coords = mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules) |> extractcoordinates
    scatter(loc1coords[1,:], loc1coords[2,:], marker=(1, stroke(0), :magenta))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(loc2coords[1,:], loc2coords[2,:], marker=(1, stroke(0), :green))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0), framestyle=:box)
end

function localizationsinsetplot_forprint(result::Result, xlims, ylims; color1 = :magenta, color2 = :green)
    loc1coords = mapreduce(x -> x.group.localizations, vcat, result.channels[1].molecules) |> extractcoordinates
    loc2coords = mapreduce(x -> x.group.localizations, vcat, result.channels[2].molecules) |> extractcoordinates
    scatter(loc1coords[1,:], loc1coords[2,:], marker=(8, stroke(0), 0.75, color1), framestyle=:none, size=(2048,2048))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(loc2coords[1,:], loc2coords[2,:], marker=(8, stroke(0), 0.75, color2))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0))
end

insetplot(result::Result, args...) = insetplot(result.channels[1].molecules, result.channels[2].molecules, args...)
function insetplot(molecules1::Vector{Molecule}, molecules2::Vector{Molecule}, xlims, ylims, include_scalebar = false)
    mol1coords = molecules1 |> extractcoordinates
    mol2coords = molecules2 |> extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(8, stroke(2, :magenta), :white))
    include_scalebar && plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(mol2coords[1,:], mol2coords[2,:], marker=(8, stroke(2, :green), :white))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0))

    loc1coords = mapreduce(x -> x.group.localizations, vcat, molecules1) |> extractcoordinates
    loc2coords = mapreduce(x -> x.group.localizations, vcat, molecules2) |> extractcoordinates
    scatter!(loc1coords[1,:], loc1coords[2,:], marker=(1, stroke(0), :magenta))
    plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
          annotations=(first(xlims) + 100, first(ylims) + 125,text("500 nm", 10, :left, :top)))

    scatter!(loc2coords[1,:], loc2coords[2,:], marker=(1, stroke(0), :green))
    plot!(aspect_ratio=:equal, xlims=xlims, yaxis=(ylims, :flip), legend=:none, grid=:hide, ticks=(0), framestyle=:box)
end

function insetplot_forprint(result::Result, xlims, ylims; color1 = :magenta, color2 = :green, include_scalebar = false)
    mol1coords = result.channels[1].molecules |> extractcoordinates
    mol2coords = result.channels[2].molecules |> extractcoordinates
    scatter(mol1coords[1,:], mol1coords[2,:], marker=(24, stroke(2, color1), 0.75, :white),
            framestyle=:none, size=(2048,2048))
    include_scalebar && plot!(first(xlims) .+ [100; 600], first(ylims) .+ [100; 100], line=(3, :black),
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

    scatter(loc1frames, loc1coords[2,:], marker=(2, :magenta, stroke(0)), xaxis=("frames", (0,22000)),
            yaxis=("y", (0,40960), :flip, 90), legend=:none)
    scatter!(loc2frames, loc2coords[2,:], marker=(2, :green, stroke(0)))
    plot!(tick_direction=:out, size=(1024,512))
end
