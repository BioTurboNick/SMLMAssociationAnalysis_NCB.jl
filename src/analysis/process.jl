
"""
    getlocalizations()

Extract `Localization`s from the data file for the given channel.
"""
function getlocalizations(alllocalizations::Vector{Localization}, channelname, startframe, nframes,
    starttrimframes, endtrimframes, radius, t_off, mergeradius)

    lowerlimit = startframe + starttrimframes
    upperlimit = startframe + nframes - endtrimframes - 1

    localizations = filter(l -> l.channel == channelname, alllocalizations)
end

"""
    getmolecules()

Group `Localization`s into `Molecule`s, using temporal and spatial proximity.
"""
function getmolecules(alllocalizations::Vector{Localization}, channelname, startframe, nframes,
    starttrimframes, endtrimframes, radius, t_off, mergeradius)

    lowerlimit = startframe + starttrimframes
    upperlimit = startframe + nframes - endtrimframes - 1

    localizations = filter(l -> l.channel == channelname, alllocalizations)
    localizationstrimmed = filter(l -> inrange(l.frame, lowerlimit, upperlimit), localizations)
    molecules = groupby_localmax_temporallimit(localizationstrimmed, radius, t_off)
    moleculesmerged = merge_close_molecules(molecules, mergeradius)

    moleculesmerged, localizationstrimmed
end
