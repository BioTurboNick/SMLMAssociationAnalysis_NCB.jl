SMLMAssociationAnalysis_NCB.jl
==============================

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4542449.svg)](https://doi.org/10.5281/zenodo.4542449)

Paper: TBD

Dataset: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4542454.svg)](https://doi.org/10.5281/zenodo.4542454)


Purpose:
--------
To examine two-channel single molecule localization microscopy (SMLM) data and quantify associations (binding)
between the two molecules.


Recommendations
---------------
Use an editor like Visual Studio Code (https://code.visualstudio.com/) for best results.


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
Follow the instructions above but use the scripts in `src/paper/`. You may obtain the data from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4542454.svg)](https://doi.org/10.5281/zenodo.4542454) and extract it to `dataset/` in the project folder.
The output files will be deposited into `output/`.
 - `run_original.jl`: read the data files, compute associations, and save result
 - `run_pos_control.jl`: read the positive control data files, compute associations, and save result
 - Analysis may be found in Pluto.jl notebooks starting with `analysis_original_`. To install Pluto, type `]add Pluto` and then `using Pluto; Pluto.run()`.
 - `run_optimize.jl`: run the calculations with various parameter choices for optimization
 - `run_optimize_analysis.jl`: generate graphs of the parameters for optimization.
 - `simulate.jl`: Simulate cells to test the algorithm.
