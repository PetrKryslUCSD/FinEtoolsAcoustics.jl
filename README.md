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

- 02/27/2024: Fix tests to pass with FinEtools 8.0.23.
- 02/25/2024: Update tutorials and documentation.
- 02/21/2024: Update for FinEtools 8.

[Past news](#past-news)

## Tutorials

There are a number of tutorials explaining the use of this package.
Check out the [index](https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl/blob/main/tutorials/index.md). 
The are tutorials themselves can be executed as
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
- Execute the desired tutorial. Here `name_tut.jl` is the name of the tutorial file:
```
include("name_tut.jl")
```

## Examples

The examples have their own environment. Change the folder to `examples`.
Then activate and instantiate the `examples` environment.
```
(FinEtoolsHeatDiff) pkg>

shell> cd examples
C:\Users\...\FinEtoolsHeatDiff.jl\examples

julia> using Pkg

julia> Pkg.activate("."); Pkg.instantiate()
  Activating project at `C:\Users\...\FinEtoolsHeatDiff.jl\examples`
   [Output suppressed...]

julia>
```

There are a number of examples, which may
be executed as described in the  [conceptual guide to
`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest).
For instance
```
julia> include("steady_state/2-d\\Poisson_examples.jl"); Poisson_examples.allrun()  
```


## <a name="past-news"></a>Past news

- 12/31/2023: Updated for Julia 1.10.0. 
- 12/17/2023: Added damping matrix for the Robin condition.
- 12/16/2023: Merge tutorials back into the package tree.
- 12/10/2023: Format source; add an example from a Fahy textbook.
- 12/08/2023: Unify terminology with the acoustics literature.
- 06/22/2023: With the exception of the transient examples, examples work.
- 06/21/2023: Updated for FinEtools 7.0.
- 05/12/2023: Updated for Julia 1.9.0. 
- 04/22/2023: Updated for generic FinEtools.
- 01/04/2023: Restructured examples. 
- 08/23/2020: Added a separate tutorial package, [FinEtoolsAcousticsTutorials.jl](https://petrkryslucsd.github.io/FinEtoolsAcousticsTutorials.jl)).
- 08/16/2020: Added tutorials.
- 08/16/2020: Dependencies updated.
- 01/23/2020: Dependencies have been updated to work with Julia 1.3.1.
- 06/11/2019: Applications are now separated  out from the `FinEtools` package.
