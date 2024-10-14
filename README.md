# Trixi simulations

This directory contains code and instructions to reproduce the numerical
experiments reported in the article

> Mimetic Metrics for the DGSEM

The results were obtained using Julia v1.10.0 and [this version of Trixi.jl](https://github.com/amrueda/Trixi.jl/tree/95603964561b2d1f9cfa83fcbd1f03935e7cc456). When this reproducibility repository was last updated, the implementation of the mimetic metrics was not yet merged into the main [Trixi.jl](https://github.com/trixi-framework/Trixi.jl/) repository.

## Instructions

To run the examples, follow the instructions:

* Move to this directory and clone the fork of Trixi.jl repository:
  ```bash
  git clone git@github.com:amrueda/Trixi.jl.git
  ```
* Move to the Trixi.jl folder and change to the branch with which the examples were computed:
  ```bash
  cd Trixi.jl
  git checkout 95603964561b2d1f9cfa83fcbd1f03935e7cc456
  ```
* Create a run direcrtoty and install all the dependencies of Trixi.jl:
  ```bash
  mkdir run
  cd run
  julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path=".."))' # Install local Trixi.jl clone
  julia --project=. -e 'using Pkg; Pkg.add(["OrdinaryDiffEq", "Trixi2Vtk", "Plots", "StaticArrays"])' # Install additional packages
  ```
* Run the examples using Julia:
  ```bash
  julia --project=. --check-bounds=no --threads=1 -e 'include(joinpath("..", "..", "tests", "elixir_advection_free_stream_mimetic_metrics.jl"))'
  ```
