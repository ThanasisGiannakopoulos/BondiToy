# BondiToy

*Toy models for Bondi-like gauges of GR. Supplementary material for
 <arxiv xxxx.xxxxx>.*


## Installation

Change to your local directory where the repository "BondiToy" is
saved. The module can be installed using Julia's REPL mode:
```
julia> ]
pkg> dev .
```

## Project Status

The module is tested against Julia 1.3.1 on Linux (Ubuntu 18.04).

## Instructions

The different toy models are under the "examples" directory labelled
as "HYPERBOLICITY_GIVENDATA_SOURCES.jl", where HYPERBOLICITY is "WH"
or "SH" for weakly and strongly hyperbolic respectively, GIVENDATA is
"smooth" or "noise" and the different SOURCES choices are "B0", "B1",
"B2" and "B3". To run the example of choice type in bash
```
julia example.jl
```

Multi-threading is supported. To enable it type in bash
```
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=n
julia example.jl
```
where n is the number of threads in multi-threading.

Data are saved in "examples/run00" and their analysis can be performed
with the scripts in "examples/data_analysis". To analyse the L2-norm
of smooth data you first need to run
"examples/data_analysis/norms_self_2d.jl" for the different
resolutions that the example with smooth data is simulated. The
resulting L2-norms can be found under
"examples/run00/example/norms_self_2".

