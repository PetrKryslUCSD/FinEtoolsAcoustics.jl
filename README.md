[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.com/PetrKryslUCSD/FinEtoolsAcoustics.jl.svg?branch=master)](https://travis-ci.com/PetrKryslUCSD/FinEtoolsAcoustics.jl)

[![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://petrkryslucsd.github.io/FinEtoolsAcoustics.jl/latest/

# FinEtoolsAcoustics: Linear acoustics application

`FinEtools` is a package for basic operations on finite element meshes.
`FinEtoolsAcoustics` is a package using `FinEtools` to solve linear acoustics problems.
Included is modal analysis, steady-state, and transient acoustics.

## News

- 06/11/2019: Applications are now separated  out from the `FinEtools` package.

[Past news](oldnews.md)

## How to run

The [FinEtools](https://github.com/PetrKryslUCSD/FinEtools.jl) package is
needed. The entire setup of `FinEtoolsAcoustics` can be performed with
```julia
] activate .; instantiate
```

The package `FinEtoolsAcoustics` can be tested as
```julia
] activate .; instantiate; test
```

There are a number of examples covering modal analysis, steady-state, and
transient acoustics. The examples may be executed as described in the  [conceptual guide to `FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest).
