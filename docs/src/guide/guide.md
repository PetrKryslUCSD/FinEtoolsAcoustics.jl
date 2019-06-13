# Guide

## Modules

The FinEtools package consists of many modules which fall into several  categories. The top-level module, `FinEtools`, includes all other modules and exports functions to constitute the public interface. The user is free to generate their own public interface, however. More details are provided [Make up your own public interface](@ref).


- Top-level:
     `FinEtoolsAcoustics` is the  top-level module.  

- Acoustics: `AlgoAcoustModule` (algorithms), `FEMMAcoustModule`, `FEMMAcoustSurfModule` (FEM machines to evaluate the matrix and vector quantities),  `MatAcoustFluidModule` (acoustic fluid material).


### Acoustics FEM machines

There is one for  the interior integrals  and one for  boundary integrals.
The  machine for the interior integrals can be used to compute:

- Evaluate the acoustic-mass matrix and the acoustic-stiffness matrix.

- Evaluate the load vector corresponding to prescribed pressure  or the prescribed second order  rate of the pressure.

The machine for the boundary integrals can be used to compute:

- Compute  transformation matrix to convert  pressure  to resultant force  or pressure to resultant torque.

- Compute the acoustic  ABC  (absorbing boundary condition) matrix.

## Algorithms

### Acoustics algorithms

At the moment there is one algorithm, for steady-state (harmonic) acoustics.

#### Example:  baffled piston

After the mesh  has been generated, the `modeldata` can be set up: Here we begin with  the region.

```julia
material = MatAcoustFluid(bulk, rho)
region1 =  FDataDict("femm"=>FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material))
```

We set up a definition of the absorbing boundary condition:

```julia
abc1  =  FDataDict("femm"=>FEMMAcoustSurf(IntegDomain(outer_fes, GaussRule(2, 2)),
          material))
```

The  surface of the piston is associated with a known-flux  boundary condition:

```julia
flux1  =  FDataDict("femm"=>FEMMAcoustSurf(IntegDomain(piston_fes, GaussRule(2, 2)),
          material),  "normal_flux"=> -rho*a_piston+0.0im);
```

And finally we make the model data,

```julia
modeldata =  FDataDict("fens"=>  fens,
                 "omega"=>omega,
                 "regions"=>[region1],
                 "flux_bcs"=>[flux1], "ABCs"=>[abc1])
```

and call  the solver:

```julia
modeldata = FinEtools.AlgoAcoustModule.steadystate(modeldata)
```

When  the algorithm completes, `modeldata["P"]` is the computed pressure field.
