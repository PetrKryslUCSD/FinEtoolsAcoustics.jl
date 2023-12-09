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

- 12/08/2023: Unify terminology with the acoustics literature.
- 06/22/2023: With the exception of the transient examples, examples work.
- 06/21/2023: Update for FinEtools 7.0.
- 05/12/2023: Updated for Julia 1.9.0. 
- 04/22/2023: Updated for generic FinEtools.
- 01/04/2023: Restructured examples. 
- 08/23/2020: Added a separate tutorial package, [FinEtoolsAcousticsTutorials.jl](https://petrkryslucsd.github.io/FinEtoolsAcousticsTutorials.jl)).
- 08/16/2020: Added tutorials.
- 08/16/2020: Dependencies updated.
- 01/23/2020: Dependencies have been updated to work with Julia 1.3.1.


[Past news](oldnews.md)

## How to test the package

Here is a record of a session to install this package and test it. You should
see something similar. The git bash running on Windows 10 was used in this
example.

Clone the repo:
```
$ git clone https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl.git
Cloning into 'FinEtoolsAcoustics.jl'...
remote: Enumerating objects: 116, done.
remote: Counting objects: 100% (116/116), done.
remote: Compressing objects: 100% (80/80), done.
remote: Total 116 (delta 37), reused 102 (delta 26), pack-reused 0
Receiving objects: 100% (116/116), 84.24 KiB | 1002.00 KiB/s, done.
Resolving deltas: 100% (37/37), done.
```
Change your working directory, and run Julia:
```
$ cd FinEtoolsAcoustics.jl/

PetrKrysl@Spectre MINGW64 /tmp/exp/FinEtoolsAcoustics.jl (master)
$ ~/AppData/Local/Julia-1.2.0-rc1/bin/julia.exe
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.2.0-rc1.0 (2019-05-30)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

```
Activate and instantiate the environment:
```
(v1.2) pkg> activate .; instantiate
[ Info: activating environment at `C:\Users\PETRKR~1\AppData\Local\Temp\exp\FinEtoolsAcoustics.jl\Project.toml`.
   Cloning default registries into `C:\Users\PetrKrysl\.julia`
   Cloning registry from "https://github.com/JuliaRegistries/General.git"
     Added registry `General` to `C:\Users\PetrKrysl\.julia\registries\General`
   Cloning git-repo `https://github.com/PetrKryslUCSD/FinEtools.jl.git`
  Updating git-repo `https://github.com/PetrKryslUCSD/FinEtools.jl.git`
 Installed SortingAlgorithms ── v0.3.1
 Installed OrderedCollections ─ v1.1.0
 Installed Arpack ───────────── v0.3.1
 Installed BinaryProvider ───── v0.5.4
 Installed DataStructures ───── v0.15.0
 Installed StatsBase ────────── v0.30.0
 Installed Missings ─────────── v0.4.1
  Building Arpack → `C:\Users\PetrKrysl\.julia\packages\Arpack\cu5By\deps\build.log`
```
Test the package:
```
(FinEtoolsAcoustics) pkg> test
   Testing FinEtoolsAcoustics
 Resolving package versions...
Test Summary: | Pass  Total
Acoustics     |   25     25
 53.500866 seconds (104.17 M allocations: 8.592 GiB, 4.88% gc time)
   Testing FinEtoolsAcoustics tests passed
```

## Examples

There are a number of examples covering modal analysis, steady-state, and
transient acoustics. The examples may be executed as described in the
[conceptual guide to `FinEtools`]
(https://petrkryslucsd.github.io/FinEtools.jl/latest).
