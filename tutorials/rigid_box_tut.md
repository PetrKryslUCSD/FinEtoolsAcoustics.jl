# Modal analysis of acoustic medium in a rigid box

Source code: [`rigid_box_tut.jl`](rigid_box_tut.jl)

## Description

Example from Boundary element acoustics: Fundamentals and computer codes, TW
Wu, page 123. Internal resonance problem. Reference frequencies: 90.7895,
181.579, 215.625, 233.959, 272.368, 281.895: computed with a quadrilateral mesh.

## Goals

- Construct the discrete model and solve the eigenvalue problem for the
  natural frequencies.
- Visualize eigenmodes.

````julia
#
````

## Definitions

The finite element code realize on the basic functionality implemented in this
package.

````julia
using FinEtools
````

The acoustics functionality is brought in:

````julia
using FinEtoolsAcoustics
````

In order to solve the eigenvalue problem we need the Arnoldi package.

````julia
import Arpack: eigs
````

The input quantities are provided without units, but the units are consistent:

````julia
rho = 1.21*1e-9;# mass density
c  = 345.0*1000;# millimeters per second
bulk =  c^2*rho;
Lx = 1900.0;# length of the box, millimeters
Ly = 800.0; # length of the box, millimeters
n = 20;#
neigvs = 18;
OmegaShift = 10.0;
````

The mesh covers a rectangle.

````julia
fens,fes = Q4block(Lx,Ly,n,n); # Mesh
````

Construct the geometry and the pressure field.

````julia
geom = NodalField(fens.xyz)
P = NodalField(zeros(size(fens.xyz,1),1))
````

There are no boundary conditions explicitly needed: The box has hard sound
boundary conditions all around. Number the unknowns in the pressure field.

````julia
numberdofs!(P)
````

The finite element machine for acoustics is constructed. The integration
rule is appropriate to the four node quadrilaterals.

````julia
femm = FEMMAcoust(IntegDomain(fes, GaussRule(2, 2)), MatAcoustFluid(bulk,rho))
````

Compute the stiffness and mass matrices:

````julia
Ka = acousticstiffness(femm, geom, P);
Ma = acousticmass(femm, geom, P);
````

And now solves eigenvalue problem. Note that the acoustic stiffness is singular,
and we need to use shifting to ensure invertibility of the first matrix.

````julia
evals, evecs, nconv = eigs(Ka+OmegaShift*Ma, Ma; nev=neigvs, which=:SM)
evals = evals .- OmegaShift;
fs = real(sqrt.(complex(evals)))./(2*pi)
println("Eigenvalues: $fs [Hz]")
````

Postprocessing: export all pressure modes.

````julia
File =  "rigid_box.vtk"
scalarllist = Any[]
for n  = 1:15
    push!(scalarllist, ("Pressure_mode_$n", deepcopy(real.(evecs[:,n]))));
end
vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.VTK.Q4;
scalars=scalarllist)
@async run(`"paraview.exe" $File`)

true
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

