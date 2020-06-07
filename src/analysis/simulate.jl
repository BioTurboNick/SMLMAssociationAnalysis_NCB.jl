"""
    simulate0()

Obtain the percentile ranks that would be expected for the set of molecules if none were bound together. This state is
simulated by randomly displacing all molecules within a given radius. This is done to simulate freedom of motion but
to roughly maintain local density. This radius should be the same as that used for montecarloaffinity().
"""
function simulate0(molecules1, molecules2, localradius, mc_iterations)
    displacements1 = randomcoordinates2d(length(molecules1), localradius)
    displacements2 = randomcoordinates2d(length(molecules2), localradius)

    newcoordinates1 = extractcoordinates(molecules1) .+ displacements1
    newcoordinates2 = extractcoordinates(molecules2) .+ displacements2

    newmolecules1 = [Molecule(Localization(i, "", newcoordinates1[1,i], newcoordinates1[2,i], newcoordinates1[3,i], 0, 1, 1)) for i ∈ 1:length(molecules1)]
    newmolecules2 = [Molecule(Localization(i, "", newcoordinates2[1,i], newcoordinates2[2,i], newcoordinates2[3,i], 0, 1, 1)) for i ∈ 1:length(molecules2)]

    neighbors1, neighbors2, distances = exclusivenearestneighbors(newmolecules1, newmolecules2)

    percentileranks = montecarloaffinity(newmolecules1, newmolecules2, neighbors1, neighbors2, distances, localradius, mc_iterations)
    percentileranks, distances
end

"""
    simulate100()

Obtain the percentile ranks that would be expected for the set of molecules of all were bound together. This state is
simulated by removing one random member of each pair and regenerating them within a given radius
of their original partner. This radius should be your best guess at how far apart fluorophores of a bound pair of your
molecule pairs should be, in nanometers. 80 nm would be a reasonable guess for average proteins with an antibody stack.
"""
function simulate100(molecules1, molecules2, neighbors1, neighbors2, bindingradius, localradius, iterations)
    replace1 = rand([true false], length(neighbors1))
    replace2 = .!replace1

    neighborindexes1_replace = [n.index for n ∈ neighbors1[replace1]]
    neighborindexes1_keep = [n.index for n ∈ neighbors1[replace2]]
    neighborindexes2_replace = [n.index for n ∈ neighbors2[replace2]]
    neighborindexes2_keep = [n.index for n ∈ neighbors2[replace1]]

    newmolecules1 = filter(m -> m.index ∉ neighborindexes1_replace, molecules1)
    newmolecules2 = filter(m -> m.index ∉ neighborindexes2_replace, molecules2)

    displacements1 = randomcoordinates2d(length(neighborindexes1_replace), bindingradius)
    displacements2 = randomcoordinates2d(length(neighborindexes2_replace), bindingradius)

    newcoordinates1 = extractcoordinates(filter(m -> m.index ∈ neighborindexes2_keep, newmolecules2)) .+ displacements1
    newcoordinates2 = extractcoordinates(filter(m -> m.index ∈ neighborindexes1_keep, newmolecules1)) .+ displacements2

    append!(newmolecules1, [Molecule(Localization(neighborindexes1_replace[i], "", newcoordinates1[1,i], newcoordinates1[2,i], newcoordinates1[3,i], 0, 1, 1)) for i ∈ 1:length(neighborindexes1_replace)])
    append!(newmolecules2, [Molecule(Localization(neighborindexes2_replace[i], "", newcoordinates2[1,i], newcoordinates2[2,i], newcoordinates2[3,i], 0, 1, 1)) for i ∈ 1:length(neighborindexes2_replace)])

    newneighbors1, newneighbors2, distances = exclusivenearestneighbors(newmolecules1, newmolecules2)

    percentileranks = montecarloaffinity(newmolecules1, newmolecules2, newneighbors1, newneighbors2, distances, localradius, iterations)
    percentileranks, distances
end
