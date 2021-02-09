"""
    Molecule

Multiple [`Localization`](@ref) objects combined into a single entity.
"""
mutable struct Molecule <: DataEntity
    index::Int
    group::LocalizationGroup
    accuracy::Float64
    Molecule(localization::Localization) =
        new(localization.index, LocalizationGroup(localization), localization.accuracy)
end

struct ResultMolecule <: DataEntity
    index::Int
    x::Float64
    y::Float64
    z::Float64
    accuracy::Float64

    ResultMolecule(mol::Molecule) = new(mol.index, mol.group.x, mol.group.y, mol.group.z, mol.accuracy)
end

import LocalizationMicroscopy.extractcoordinates
extractcoordinates(molecule::Molecule) = [molecule.group.x; molecule.group.y; molecule.group.z]

import Base.push!
function push!(molecule::Molecule, newlocalization::Localization...)
    push!(molecule.group, newlocalization...)
    calcaccuracy!(molecule)
end

function push!(molecule::Molecule, othermolecule::Molecule)
    append!(molecule.group, othermolecule.group)
end

import Base.append!
function append!(molecule::Molecule, newlocalizations::Vector{Localization})
    append!(molecule.group, newlocalizations)
    calcaccuracy!(molecule)
end

function calcaccuracy!(molecule::Molecule)
    count = length(molecule.group.localizations)
    molecule.accuracy = mean(l -> l.accuracy, molecule.group.localizations) / sqrt(count)
end

export Molecule, ResultMolecule, push!, append!
