# Trixi simulations

This directory contains code and instructions to reproduce the numerical
experiments reported in the article

* Bach, D., Rueda-Ramírez, A., Kopriva, D. A., & Gassner, G. J. (2024). Mimetic Metrics for the DGSEM. [arXiv preprint:2410.14502](https://arxiv.org/abs/2410.14502).

The results were obtained using Julia v1.11.5 and [this version of Trixi.jl](https://github.com/amrueda/Trixi.jl/tree/987cdbb78924d1159dc14d2e6e78de4a31e94770) (release tag [mimetic_metrics_paper](https://github.com/amrueda/Trixi.jl/releases/tag/mimetic_metrics_paper)). When this reproducibility repository was last updated, the implementation of the mimetic metrics was not yet merged into the main [Trixi.jl](https://github.com/trixi-framework/Trixi.jl/) repository.

## Instructions

To run the examples, follow the instructions:

* Clone the reproducibility repository and move to the cloned directory:
   ```bash
  git clone git@github.com:amrueda/paper_2024_mimetic_metrics.git
  cd paper_2024_mimetic_metrics
  ```
* Clone our fork of Trixi.jl, move to the Trixi.jl folder, and change to the branch with which the examples were computed (you can alternatively download the [tagged source code](https://github.com/amrueda/Trixi.jl/releases/tag/mimetic_metrics_paper)):
  ```bash
  git clone git@github.com:amrueda/Trixi.jl.git
  cd Trixi.jl
  git checkout 987cdbb78924d1159dc14d2e6e78de4a31e94770
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
  julia --project=. --threads=1 -e 'include(joinpath("..", "..", "tests", "elixir_advection_free_stream_mimetic_metrics.jl"))'
  julia --project=. --threads=1 -e 'include(joinpath("..", "..", "tests", "elixir_euler_free_stream_mimetic_metrics.jl"))'
  ```
* We also have a light-weight example (not shown in the paper) that runs the linear advection equation on a curvilinear grid
  ```bash
  julia --project=. --threads=1 -e 'include(joinpath("..", "..", "tests", "elixir_advection_free_stream_mimetic_metrics.jl"))'
  julia --project=. --threads=1 -e 'include(joinpath("..", "..", "tests", "elixir_euler_free_stream_mimetic_metrics.jl"))'
  ```
