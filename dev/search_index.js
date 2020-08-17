var documenterSearchIndex = {"docs":
[{"location":"tutorials/rigid_box_tut.html#Modal-analysis-of-acoustic-medium-in-a-rigid-box","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"","category":"section"},{"location":"tutorials/rigid_box_tut.html#Description","page":"Modal analysis of acoustic medium in a rigid box","title":"Description","text":"","category":"section"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 123. Internal resonance problem. Reference frequencies: 90.7895, 181.579, 215.625, 233.959, 272.368, 281.895 Quadrilateral mesh.","category":"page"},{"location":"tutorials/rigid_box_tut.html#Goals","page":"Modal analysis of acoustic medium in a rigid box","title":"Goals","text":"","category":"section"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"Introduce definition of cross-section.\nShow generation of finite element mesh of beams.\nDescribe geometry, displacement, and rotation fields.\nDescribe application of support conditions.\nCalculate the discrete model quantities and solve the free vibration problem.\nDemonstrate visualization of the free vibrations.","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"#","category":"page"},{"location":"tutorials/rigid_box_tut.html#Definitions","page":"Modal analysis of acoustic medium in a rigid box","title":"Definitions","text":"","category":"section"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"The finite element code realize on the basic functionality implemented in this package.","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"using FinEtools","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"The acoustics functionality is brought in:","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"using FinEtoolsAcoustics","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"In order to solve the eigenvalue problem we need the Arnoldi package.","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"import Arpack: eigs","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"The input quantities are provided without units, but the units are consistent:","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"rho = 1.21*1e-9;# mass density\nc  = 345.0*1000;# millimeters per second\nbulk =  c^2*rho;\nLx = 1900.0;# length of the box, millimeters\nLy = 800.0; # length of the box, millimeters\nn = 14;#\nneigvs = 18;\nOmegaShift = 10.0;","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"The mesh covers a rectangle.","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"fens,fes = Q4block(Lx,Ly,n,n); # Mesh","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"Construct the geometry and the pressure field.","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"geom = NodalField(fens.xyz)\nP = NodalField(zeros(size(fens.xyz,1),1))","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"There are no boundary conditions explicitly needed: The box has hard sound boundary conditions all around. Number the unknowns in the pressure field.","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"numberdofs!(P)","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"The finite element machine for acoustics is constructed. The integration rule is appropriate to the four node quadrilaterals.","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"femm = FEMMAcoust(IntegDomain(fes, GaussRule(2, 2)), MatAcoustFluid(bulk,rho))","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"Compute the stiffness and mass matrices:","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"S = acousticstiffness(femm, geom, P);\nC = acousticmass(femm, geom, P);","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"And now solves eigenvalue problem. Note that the acoustic mass is singular, and we need to use shifting to ensure convertibility of the first matrix.","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"evals, evecs, nconv = eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)\nevals = evals .- OmegaShift;\nfs = real(sqrt.(complex(evals)))./(2*pi)\nprintln(\"Eigenvalues: $fs [Hz]\")","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"Postprocessing: export all pressure modes.","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"File =  \"rigid_box.vtk\"\nscalarllist = Any[]\nfor n  = 1:15\n    push!(scalarllist, (\"Pressure_mode_$n\", evecs[:,n]));\nend\nvtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.VTK.Q4;\nscalars=scalarllist)\n@async run(`\"paraview.exe\" $File`)\n\ntrue","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"","category":"page"},{"location":"tutorials/rigid_box_tut.html","page":"Modal analysis of acoustic medium in a rigid box","title":"Modal analysis of acoustic medium in a rigid box","text":"This page was generated using Literate.jl.","category":"page"},{"location":"guide/guide.html#Guide","page":"How to guide","title":"Guide","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The FinEtools package is used here to solve heat conduction problems.","category":"page"},{"location":"guide/guide.html#Modules","page":"How to guide","title":"Modules","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The package FinEtoolsAcoustics has the following structure:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"Top-level:     FinEtoolsAcoustics is the  top-level module.  \nAcoustics: AlgoAcoustModule (algorithms), FEMMAcoustModule, FEMMAcoustSurfModule (FEM machines to evaluate the matrix and vector quantities),  MatAcoustFluidModule (acoustic fluid material).","category":"page"},{"location":"guide/guide.html#Acoustics-FEM-machines","page":"How to guide","title":"Acoustics FEM machines","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"There is one for  the interior integrals  and one for  boundary integrals. The  machine for the interior integrals can be used to compute:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"Evaluate the acoustic-mass matrix and the acoustic-stiffness matrix.\nEvaluate the load vector corresponding to prescribed pressure  or the prescribed second order  rate of the pressure.","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The machine for the boundary integrals can be used to compute:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"Compute  transformation matrix to convert  pressure  to resultant force  or pressure to resultant torque.\nCompute the acoustic  ABC  (absorbing boundary condition) matrix.","category":"page"},{"location":"guide/guide.html#Acoustics-algorithms","page":"How to guide","title":"Acoustics algorithms","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"At the moment there is one algorithm, for steady-state (harmonic) acoustics.","category":"page"},{"location":"guide/guide.html#Example:-baffled-piston","page":"How to guide","title":"Example:  baffled piston","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"After the mesh  has been generated, the modeldata can be set up: Here we begin with  the region.","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"material = MatAcoustFluid(bulk, rho)\nregion1 =  FDataDict(\"femm\"=>FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material))","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"We set up a definition of the absorbing boundary condition:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"abc1  =  FDataDict(\"femm\"=>FEMMAcoustSurf(IntegDomain(outer_fes, GaussRule(2, 2)),\n          material))","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The  surface of the piston is associated with a known-flux  boundary condition:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"flux1  =  FDataDict(\"femm\"=>FEMMAcoustSurf(IntegDomain(piston_fes, GaussRule(2, 2)),\n          material),  \"normal_flux\"=> -rho*a_piston+0.0im);","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"And finally we make the model data,","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"modeldata =  FDataDict(\"fens\"=>  fens,\n                 \"omega\"=>omega,\n                 \"regions\"=>[region1],\n                 \"flux_bcs\"=>[flux1], \"ABCs\"=>[abc1])","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"and call  the solver:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"modeldata = FinEtools.AlgoAcoustModule.steadystate(modeldata)","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"When  the algorithm completes, modeldata[\"P\"] is the computed pressure field.","category":"page"},{"location":"tutorials/baffled_piston_tut.html#Baffled-piston-in-a-half-sphere-domain","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"","category":"section"},{"location":"tutorials/baffled_piston_tut.html#Description","page":"Baffled piston in a half-sphere domain","title":"Description","text":"","category":"section"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Baffled piston in a half-sphere domain with Absorbing Boundary Condition (ABC).","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"(Image: )","category":"page"},{"location":"tutorials/baffled_piston_tut.html#Goals","page":"Baffled piston in a half-sphere domain","title":"Goals","text":"","category":"section"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Demonstrate the use of an algorithm to run the simulation.","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"#","category":"page"},{"location":"tutorials/baffled_piston_tut.html#Definitions","page":"Baffled piston in a half-sphere domain","title":"Definitions","text":"","category":"section"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"The finite element code relies on the basic functionality implemented in this package.","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"using FinEtools","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"The acoustics functionality is brought in:","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"using FinEtoolsAcoustics","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"We shall need some facilities from the linear algebra package","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"using LinearAlgebra","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"The input quantities are provided including the units:","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"rho = 1.21*phun(\"kg/m^3\");# mass density\nc  = 343.0*phun(\"m/s\");# sound speed\nbulk =  c^2*rho;\nomega =  7500*phun(\"rev/s\");      # frequency of the piston\na_piston =  -1.0*phun(\"mm/s\")     # amplitude of the piston acceleration\nR = 50.0*phun(\"mm\");# radius of the piston\nRo = 150.0*phun(\"mm\"); # radius of the enclosure\nnref = 4;#number of refinements of the sphere around the piston\nnlayers = 35;                     # number of layers of elements surrounding the piston\ntolerance = R/(2^nref)/100","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Hexahedral mesh","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"fens,fes  =  H8sphere(R,nref);\nbfes  =  meshboundary(fes)\nFile  =   \"baffledabc_boundary.vtk\"\nvtkexportmesh(File, connasarray(bfes), fens.xyz, FinEtools.MeshExportModule.VTK.Q4)\n@async run(`\"paraview.exe\" $File`)\n\nl = selectelem(fens, bfes, facing = true, direction = [1.0 1.0  1.0], dotmin= 0.001)\nex(xyz, layer) = (R+layer/nlayers*(Ro-R))*xyz/norm(xyz)\nfens1, fes1  =  H8extrudeQ4(fens, subset(bfes,l), nlayers, ex);\nfens, newfes1, fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)\nfes = cat(newfes1,fes2)","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Piston surface mesh","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"bfes  =  meshboundary(fes)\nl1 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])\nl2 = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0], inflate = tolerance)\npiston_fes = subset(bfes,intersect(l1,l2));","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Outer spherical boundary","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"louter = selectelem(fens, bfes, facing = true, direction = [1.0 1.0  1.0], dotmin= 0.001)\nouter_fes = subset(bfes,louter);","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"The simulation is driven by setting up and then executing an algorithm. For this purpose we store the data in a \"model data dictionary\".","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Region of the fluid. Define the finite element machine for the fluid region.","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"material = MatAcoustFluid(bulk, rho)\nregion1 =  FDataDict(\"femm\"=>FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material))","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Define the finite element machine for the ABC (absorbing boundary condition).","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"abc1  =  FDataDict(\"femm\"=>FEMMAcoustSurf(IntegDomain(outer_fes, GaussRule(2, 2)), material))","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"The normal flux is prescribed on the surface of the piston in terms of the known acceleration.","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"flux1  =  FDataDict(\"femm\"=>FEMMAcoustSurf(IntegDomain(piston_fes, GaussRule(2, 2)), material),  \"normal_flux\"=> -rho*a_piston+0.0im);","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Make model data dictionary. It completely defines the parameters of the problem.","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"modeldata =  FDataDict(\"fens\"=>fens, \"omega\"=>omega, \"regions\"=>[region1], \"flux_bcs\"=>[flux1], \"ABCs\"=>[abc1])","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Call the solver. The model data is returned enriched of the solution parameters, such as the pressure field P.","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"modeldata = FinEtoolsAcoustics.AlgoAcoustModule.steadystate(modeldata)","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Extract geometry field and the pressure field for postprocessing.","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"geom = modeldata[\"geom\"]\nP = modeldata[\"P\"]","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"Write out the postprocessing data is a VTK file. Dump both the magnitude and the components of the solution (complex field).","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"File  =   \"baffledabc.vtk\"\nvtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.VTK.H8; scalars = [(\"absP\", abs.(P.values)), (\"realP\", real.(P.values)), (\"imagP\", imag.(P.values))])","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"If the paraview program is installed, run it.","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"@async run(`\"paraview.exe\" $File`)\n\ntrue","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"","category":"page"},{"location":"tutorials/baffled_piston_tut.html","page":"Baffled piston in a half-sphere domain","title":"Baffled piston in a half-sphere domain","text":"This page was generated using Literate.jl.","category":"page"},{"location":"tutorials/tutorials.html","page":"Tutorials","title":"Tutorials","text":"Table of contents","category":"page"},{"location":"tutorials/tutorials.html#Tutorials","page":"Tutorials","title":"Tutorials","text":"","category":"section"},{"location":"tutorials/tutorials.html","page":"Tutorials","title":"Tutorials","text":"Modes in a rigid box\nAcoustic field of a baffled piston","category":"page"},{"location":"index.html#FinEtoolsAcoustics-Documentation","page":"Home","title":"FinEtoolsAcoustics Documentation","text":"","category":"section"},{"location":"index.html#Tutorials","page":"Home","title":"Tutorials","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The tutorials are provided in the form of Julia scripts and Markdown files. ","category":"page"},{"location":"index.html#Manual","page":"Home","title":"Manual","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The description of the types and the functions, organized by module and/or other logical principle.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n    \"man/types.md\",\n    \"man/functions.md\",\n]\nDepth = 2","category":"page"},{"location":"index.html#Conceptual-guide","page":"Home","title":"Conceptual guide","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The construction of the toolkit is described: the composition of modules, the basic data structures, the methodology of computing quantities required in the finite element methodology, and more.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n    \"guide/guide.md\",\n]\nDepth = 1","category":"page"},{"location":"man/types.html#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"man/types.html#FEM-machines","page":"Types","title":"FEM machines","text":"","category":"section"},{"location":"man/types.html#Acoustics","page":"Types","title":"Acoustics","text":"","category":"section"},{"location":"man/types.html","page":"Types","title":"Types","text":"Modules = [FinEtools, FinEtoolsAcoustics.FEMMAcoustModule, FinEtoolsAcoustics.FEMMAcoustSurfModule]\nPrivate = true\nOrder = [:type]","category":"page"},{"location":"man/types.html#FinEtoolsAcoustics.FEMMAcoustModule.FEMMAcoust","page":"Types","title":"FinEtoolsAcoustics.FEMMAcoustModule.FEMMAcoust","text":"FEMMAcoust{S<:AbstractFESet, F<:Function, M} <: AbstractFEMM\n\nType for linear acoustics finite element modeling machine.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.FEMMAcoustSurf","page":"Types","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.FEMMAcoustSurf","text":"FEMMAcoustSurf{S<:AbstractFESet, F<:Function, M, NF<:Function} <: AbstractFEMM\n\nClass for linear acoustics finite element modeling machine.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.FEMMAcoustSurf-Union{Tuple{M}, Tuple{F}, Tuple{S}, Tuple{IntegDomain{S,F},M}} where M where F<:Function where S<:AbstractFESet","page":"Types","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.FEMMAcoustSurf","text":"FEMMAcoustSurf(integdomain::IntegDomain{S, F},  material::M) where {S<:AbstractFESet, F<:Function, M}\n\nCreate the FEMM for integrals over the surface. The normal is computed from the geometry of the surface elements.\n\n\n\n\n\n","category":"method"},{"location":"man/types.html#Material-models","page":"Types","title":"Material models","text":"","category":"section"},{"location":"man/types.html#Material-models-for-acoustics","page":"Types","title":"Material models for acoustics","text":"","category":"section"},{"location":"man/types.html","page":"Types","title":"Types","text":"Modules = [FinEtools, FinEtools.MatModule, FinEtoolsAcoustics.MatAcoustFluidModule]\nPrivate = true\nOrder = [:type]","category":"page"},{"location":"man/types.html#FinEtools.MatModule.AbstractMat","page":"Types","title":"FinEtools.MatModule.AbstractMat","text":"AbstractMat\n\nAbstract type of material.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsAcoustics.MatAcoustFluidModule.MatAcoustFluid","page":"Types","title":"FinEtoolsAcoustics.MatAcoustFluidModule.MatAcoustFluid","text":"MatAcoustFluid <: AbstractMat\n\nType for acoustic fluid material.\n\n\n\n\n\n","category":"type"},{"location":"man/functions.html#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"man/functions.html#FEM-machines","page":"Functions","title":"FEM machines","text":"","category":"section"},{"location":"man/functions.html#Acoustics","page":"Functions","title":"Acoustics","text":"","category":"section"},{"location":"man/functions.html","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsAcoustics.FEMMAcoustModule, FinEtoolsAcoustics.FEMMAcoustSurfModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustModule.acousticmass-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoust,A,NodalField,NodalField{T}}} where A<:AbstractSysmatAssembler where T<:Number","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustModule.acousticmass","text":"acousticmass(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the acoustic mass matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler\ngeom = geometry field\nP = acoustic (perturbation) pressure field\n\nReturn a matrix.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustModule.acousticstiffness-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoust,A,NodalField,NodalField{T}}} where A<:AbstractSysmatAssembler where T<:Number","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustModule.acousticstiffness","text":"acousticstiffness(self::FEMMAcoust, assembler::A,\n  geom::NodalField,\n  Pddot::NodalField{T}) where {T<:Number,\n  A<:AbstractSysmatAssembler}\n\nCompute the acoustic stiffness matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler\ngeom = geometry field\nPddot = second order rate of the acoustic (perturbation) pressure field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustModule.nzebcloadsacousticmass-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoust,A,NodalField,NodalField{T}}} where A<:AbstractSysvecAssembler where T<:Number","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustModule.nzebcloadsacousticmass","text":"nzebcloadsacousticmass(self::FEMMAcoust, assembler::A,\n  geom::NodalField, P::NodalField{T}) where {T<:Number,\n  A<:AbstractSysvecAssembler}\n\nCompute load vector for nonzero EBC for prescribed pressure.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler\ngeom = geometry field\nP = acoustic (perturbation) pressure field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustModule.nzebcloadsacousticstiffness-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoust,A,NodalField,NodalField{T}}} where A<:AbstractSysvecAssembler where T<:Number","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustModule.nzebcloadsacousticstiffness","text":"nzebcloadsacousticstiffness(self::FEMMAcoust, assembler::A,\n  geom::NodalField,\n  Pddot::NodalField{T}) where {T<:Number,\n  A<:AbstractSysvecAssembler}\n\nCompute load vector for nonzero EBC for prescribed second-order pressure rate.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler\ngeom = geometry field\nPddot = second order rate of the acoustic (perturbation) pressure field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticABC-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoustSurf,A,NodalField,NodalField{T}}} where A<:AbstractSysmatAssembler where T<:Number","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticABC","text":"acousticABC(self::FEMMAcoustSurf, assembler::A,\n  geom::NodalField,\n  Pdot::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the acoustic ABC (Absorbing Boundary Condition) matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler; must be able to assemble unsymmetric matrix\ngeom = geometry field\nPdot = rate of the acoustic (perturbation) pressure field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticcouplingpanels-Union{Tuple{T}, Tuple{A}, Tuple{FEMMAcoustSurf,A,NodalField,NodalField{T}}} where T where A<:AbstractSysmatAssembler","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticcouplingpanels","text":"acousticcouplingpanels(self::FEMMAcoustSurf, geom::NodalField, u::NodalField{T}) where {T}\n\nCompute the acoustic pressure-structure coupling matrix.\n\nThe acoustic pressure-nodal force matrix transforms the pressure distributed along the surface to forces acting on the nodes of the finite element model. Its transpose transforms displacements (or velocities, or accelerations) into the normal component of the displacement (or velocity, or acceleration) along the surface.\n\nArguments\n\ngeom=geometry field\nu = displacement field\n\nnote: Note\n\n\nn = outer normal (pointing into the acoustic medium).\nThe pressures along the surface are assumed constant (uniform) along   each finite element –- panel. The panel pressures are assumed to be   given the same numbers as the serial numbers of the finite elements in   the set.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultantforce-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoustSurf,A,NodalField,NodalField{T},GeneralField}} where A<:AbstractSysmatAssembler where T<:Number","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultantforce","text":"pressure2resultantforce(self::FEMMAcoustSurf, assembler::A,\n  geom::NodalField,\n  P::NodalField{T},\n   Force::Field) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the rectangular coupling matrix that transcribes given pressure on the surface into the resultant force acting on the surface.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler; must be able to assemble unsymmetric matrix\ngeom = geometry field\nP = acoustic (perturbation) pressure field\nForce = field for the force resultant\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultanttorque-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoustSurf,A,NodalField,NodalField{T},GeneralField,Array{Float64,1}}} where A<:AbstractSysmatAssembler where T<:Number","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultanttorque","text":"pressure2resultanttorque(self::FEMMAcoustSurf, assembler::A,\n  geom::NodalField,\n  P::NodalField{T},\n  Torque::GeneralField, CG::FFltVec) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the rectangular coupling matrix that transcribes given pressure on the surface into the resultant torque acting on the surface with respect to the CG.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler; must be able to assemble unsymmetric matrix\ngeom = geometry field\nP = acoustic (perturbation) pressure field\nTorque = field for the torque resultant\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Algorithms","page":"Functions","title":"Algorithms","text":"","category":"section"},{"location":"man/functions.html#Acoustics-2","page":"Functions","title":"Acoustics","text":"","category":"section"},{"location":"man/functions.html","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsAcoustics.AlgoAcoustModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsAcoustics.AlgoAcoustModule.steadystate-Tuple{Dict{String,Any}}","page":"Functions","title":"FinEtoolsAcoustics.AlgoAcoustModule.steadystate","text":"steadystate(modeldata::FDataDict)\n\nSteady-state acoustics solver.\n\nmodeldata = dictionary with string keys\n\n\"fens\"  = finite element node set\n\"regions\"  = array of region dictionaries\n\"essential_bcs\" = array of essential boundary condition dictionaries\n\"ABCs\" = array of absorbing boundary condition dictionaries\n\"flux_bcs\" = array of flux boundary condition dictionaries\n\nFor each region (connected piece of the domain made of a particular material), mandatory, the  region dictionary  contains items:\n\n\"femm\" = finite element mmodel machine (mandatory);\n\nFor essential boundary conditions (optional) each dictionary would hold\n\n\"pressure\" = fixed (prescribed) pressure (scalar),  or         a function with signature             function T = f(x)         If not given, zero pressure assumed.\n\"node_list\" = list of nodes on the boundary to which the condition applies         (mandatory)\n\nFor absorbing boundary conditions (optional) each dictionary may hold\n\n\"femm\" = finite element mmodel machine (mandatory).\n\nFor flux boundary conditions (optional) each dictionary would hold\n\n\"femm\" = finite element mmodel machine (mandatory);\n\"normal_flux\" = normal component of the flux through the boundary (scalar),     which is the normal derivative of the pressure.\n\nOutput\n\nmodeldata = the dictionary is augmented with\n\n\"geom\" = the nodal field that is the geometry\n\"P\" = the nodal field that is the computed pressure (in the general a           complex-number field)\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Material-models","page":"Functions","title":"Material models","text":"","category":"section"},{"location":"man/functions.html#Material-models-for-acoustics","page":"Functions","title":"Material models for acoustics","text":"","category":"section"},{"location":"man/functions.html","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsAcoustics.MatAcoustFluidModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsAcoustics.MatAcoustFluidModule.bulkmodulus-Tuple{MatAcoustFluid}","page":"Functions","title":"FinEtoolsAcoustics.MatAcoustFluidModule.bulkmodulus","text":"bulkmodulus(self::MatAcoustFluid)\n\nReturn the bulk modulus.\n\n\n\n\n\n","category":"method"}]
}
