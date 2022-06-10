This repository is a [julia](https://github.com/JuliaLang/julia) implementation of the model described in the group project "A Static Networks Approach to Exploring Collective Animal Motion".

Specifically, the information propagation within that model is computationally explored here.

# Usage
First, ensure you have `julia` installed and executable added to your PATH. Then, install the relevant packages.
Note that the package [PGFPlotsX.jl](https://github.com/kristofferc/PGFPlotsX.jl) may require active LaTeX installation.

Model parameters are specified under the directory `profile/`.


Use the following command to run the simulation.
```shell
julia fish-school-run.jl Profile_Name [overwrite]
```

To run all the profiles, use the following command.
```shell
julia fish-school-compile.jl
```


# Figures
Source files of figures used in the paper are placed under the directory `fig/`, with name `fig_result_*.tex`. They rely on the data files under the same directory and can be compiled using LuaLaTeX.

