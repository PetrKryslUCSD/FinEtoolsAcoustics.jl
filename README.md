[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.com/PetrKryslUCSD/FinEtoolsAcoustics.jl.svg?branch=master)](https://travis-ci.com/PetrKryslUCSD/FinEtoolsAcoustics.jl)

[![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://petrkryslucsd.github.io/FinEtoolsAcoustics.jl/latest/

# FinEtoolsAcoustics: Linear acoustics application

[`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl.git) is a package
for basic operations on finite element meshes. `FinEtoolsHeatDiff` is a package
using `FinEtools` to solve linear heat conduction (diffusion) problems.

## News

- 06/11/2019: Applications are now separated  out from the `FinEtools` package.

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
[conceptual guide to
`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest).
