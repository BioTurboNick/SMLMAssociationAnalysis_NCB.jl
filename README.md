SMLMAssociationAnalysis_NCB.jl
==============================

Purpose:
--------
To examine two-channel single molecule localization microscopy (SMLM) data and quantify associations (binding)
between them.

Recommendations
---------------
Use an editor like Juno (Atom) for best results.


Use (arbitrary data)
--------------------
1. Adapt `run.jl` to your own files by modifying the parameters in the top section to suit your files.
2. Load the module:
```
pkg> activate .
julia> using SMLMAssociationAnalysis_NCB
```
(`pkg>` prompt is accessed by typing the ']' first; backspace to return to `julia>` prompt.)
(Alternatively, if you are launching julia from a terminal, navigate to the project directory and execute
`julia --project`)
3. Execute `run.jl`. (If using the console only, `include("run.jl")` will execute the script.)

Recrate the paper's data
------------------------
The underlying data files have some quirks that need to be adjusted, so a special file is provided.
Download the data from `<repository>`.
Follow the instructions above but execute `run_original.jl` instead, followed by `run_plots.jl`
