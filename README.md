SMLMAssociationAnalysis_NCB.jl
==============================

Purpose:
--------
To examine two-channel single molecule localization microscopy (SMLM) data and quantify associations (binding)
between the two molecules.

Recommendations
---------------
Use an editor like Juno (https://junolab.org/; or Atom with `uber-juno` plugin) for best results. Visual Studio Code has
a Julia plugin but is missing some features.


Use for your own data
---------------------
1. Adapt `src/run.jl` to your own files by modifying the parameters in the top section to suit your files.
2. Load the module:
```
pkg> activate .                          # sets the current path as the environment to use
julia> using SMLMAssociationAnalysis_NCB # load the module
```
(`pkg>` prompt is accessed by typing ']' first; backspace to return to `julia>` prompt.)
(Alternatively, if you are launching julia from a terminal, navigate to the project directory and execute
`julia --project` instead)

3. Execute `run.jl`. (If using the console only, `include("run.jl")` will execute the script.)
4. Execute `run_plots.jl` to generate some basic plots.

Recrate the paper's analysis
----------------------------
The underlying data files have some quirks that need to be adjusted, so a special run file is provided.
Follow the instructions above but use the scripts in `src/paper/`. The data is contained in `dataset/` and will be
deposited into `output/`.
 - run_original.jl: read the data files, compute associations, and save result
 - run_original_analysis.jl: load result and calculate statistics and plots
 - run_optimize.jl: run the calculations with various parameter choices for optimization
 - run_optimize_analysis.jl: generate graphs of the parameters for optimization.
 - simulate.jl: Simulate cells to test the algorithm.

Key to dataset:
     Mdm2-p53: A = + Dox / - Nut
               B = - Dox / - Nut
               C = + Dox / + Nut
               D = - Dox / + Nut
     MEG3-p53: A = + Dox / MEG3
               B = - Dox / MEG3
               C = + Dox / GAPDH
               D = - Dox / GAPDH
