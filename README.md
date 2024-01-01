[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl/actions)
[![Code Coverage](https://codecov.io/gh/PetrKryslUCSD/FinEtoolsAcoustics.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/PetrKryslUCSD/FinEtoolsAcoustics.jl)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtoolsAcoustics.jl/latest)
[![Codebase Graph](https://img.shields.io/badge/Codebase-graph-green.svg)](https://octo-repo-visualization.vercel.app/?repo=PetrKryslUCSD/FinEtoolsAcoustics.jl)


# FinEtoolsAcoustics: Linear acoustics application

[`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl.git) is a package
for basic operations on finite element meshes. `FinEtoolsAcoustics` is a package
using `FinEtools` to solve linear acoustics problems.

## News

- 12/31/2023: Updated for Julia 1.10.0. 
- 12/17/2023: Added damping matrix for the Robin condition.
- 12/16/2023: Merge tutorials back into the package tree.

[Past news](oldnews.md)

## Tutorials

There are a number of tutorials explaining the use of this package.
Check out the [index](https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl/blob/main/tutorials/index.md). The are tutorials themselves can be executed as
follows:

- Download the package or clone it.
```
git clone https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl.git
```
- Change into the `tutorials` folder: `cd .\FinEtoolsAcoustics.jl\tutorials`.
- Start Julia: `julia`.
- Activate the environment:
```
using Pkg; Pkg.activate("."); Pkg.instantiate();
```
- Execute the desired tutorial. Here `name.jl` is the name of the tutorial file:
```
include("name.jl")
```

## Examples

There are a number of examples covering modal analysis, steady-state, and
transient acoustics. The examples may be executed as described in the
[conceptual guide to `FinEtools`]
(https://petrkryslucsd.github.io/FinEtools.jl/latest).
