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
()
3. Execute `run.jl`. (If using the console only, `include("run.jl")` will execute the script.)

Use (recreate analysis)
-----------------------
The underlying data files have some quirks that need to be adjusted, so a special file is provided.
Follow the instructions above but execute `run_original.jl` instead.
