var documenterSearchIndex = {"docs":
[{"location":"guide/guide.html#Guide","page":"How to guide","title":"Guide","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The FinEtools package is used here to solve linear acoustics problems.","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"Tutorials  are provided in the form of Julia scripts and Markdown files in a dedicated folder: index of tutorials. ","category":"page"},{"location":"guide/guide.html#Modules","page":"How to guide","title":"Modules","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The package FinEtoolsAcoustics has the following structure:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"Top-level:     FinEtoolsAcoustics is the  top-level module.  \nAcoustics: AlgoAcoustModule (algorithms), FEMMAcoustModule, FEMMAcoustSurfModule (FEM machines to evaluate the matrix and vector quantities),  MatAcoustFluidModule (acoustic fluid material).","category":"page"},{"location":"guide/guide.html#Acoustics-FEM-machines","page":"How to guide","title":"Acoustics FEM machines","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"There is one for  the interior integrals  and one for  boundary integrals. The  machine for the interior integrals can be used to compute:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"Evaluate the acoustic-mass matrix and the acoustic-stiffness matrix.","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The machine for the boundary integrals can be used to compute:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"Compute  transformation matrix to convert  pressure  to resultant force  or pressure to resultant torque.\nCompute the acoustic  ABC  (absorbing boundary condition) matrix.\nCompute the acoustic Robin  boundary condition matrix.","category":"page"},{"location":"guide/guide.html#Acoustics-algorithms","page":"How to guide","title":"Acoustics algorithms","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"At the moment there is one algorithm, for steady-state (harmonic) acoustics.","category":"page"},{"location":"guide/guide.html#Example:-baffled-piston","page":"How to guide","title":"Example:  baffled piston","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"After the mesh  has been generated, the modeldata can be set up: Here we begin with  the region.","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"material = MatAcoustFluid(bulk, rho)\nregion1 =  FDataDict(\"femm\"=>FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material))","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"We set up a definition of the absorbing boundary condition:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"abc1  =  FDataDict(\"femm\"=>FEMMAcoustSurf(IntegDomain(outer_fes, GaussRule(2, 2)),\n          material))","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The  surface of the piston is associated with a known-flux  boundary condition (prescribed normal component of the velocity):","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"flux1  =  FDataDict(\"femm\"=>FEMMAcoustSurf(IntegDomain(piston_fes, GaussRule(2, 2)),\n          material),  \"normal_flux\"=> -1.0im*rho*omega*v_piston);","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"And finally we make the model data,","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"modeldata =  FDataDict(\"fens\"=>  fens,\n                 \"omega\"=>omega,\n                 \"regions\"=>[region1],\n                 \"flux_bcs\"=>[flux1], \"ABCs\"=>[abc1])","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"and call  the solver:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"modeldata = FinEtools.AlgoAcoustModule.steadystate(modeldata)","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"When  the algorithm completes, modeldata[\"P\"] is the computed pressure field.","category":"page"},{"location":"index.html#FinEtoolsAcoustics-Documentation","page":"Home","title":"FinEtoolsAcoustics Documentation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"","category":"page"},{"location":"index.html#Conceptual-guide","page":"Home","title":"Conceptual guide","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The construction of the toolkit is described: the composition of modules, the basic data structures, the methodology of computing quantities required in the finite element methodology, and more.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n    \"guide/guide.md\",\n]\nDepth = 1","category":"page"},{"location":"index.html#Manual","page":"Home","title":"Manual","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The description of the types and the functions, organized by module and/or other logical principle.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n    \"man/man.md\",\n]\nDepth = 2","category":"page"},{"location":"index.html#Tutorials","page":"Home","title":"Tutorials","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The tutorials are provided in the form of Julia scripts and Markdown files in the tutorials folder. ","category":"page"},{"location":"index.html#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"","category":"page"},{"location":"man/man.html#Manual","page":"Manual","title":"Manual","text":"","category":"section"},{"location":"man/man.html","page":"Manual","title":"Manual","text":"CurrentModule = FinEtoolsAcoustics","category":"page"},{"location":"man/man.html#FEM-machines","page":"Manual","title":"FEM machines","text":"","category":"section"},{"location":"man/man.html#Acoustics:-volume","page":"Manual","title":"Acoustics: volume","text":"","category":"section"},{"location":"man/man.html","page":"Manual","title":"Manual","text":"FEMMAcoust\nacousticstiffness\nacousticmass\ninspectintegpoints","category":"page"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustModule.FEMMAcoust","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustModule.FEMMAcoust","text":"FEMMAcoust{ID<:IntegDomain, M} <: AbstractFEMM\n\nType for linear acoustics finite element modeling machine.\n\n\n\n\n\n","category":"type"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustModule.acousticstiffness","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustModule.acousticstiffness","text":"acousticstiffness(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the acoustic \"stiffness\" matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler\ngeom = geometry field\nP = acoustic (perturbation) pressure field\n\nReturn a matrix.\n\nThe acoustic \"stiffness\" matrix is by convention called stiffness, however its mechanical meaning is quite different. (It has to do with kinetic energy.) It is the matrix mathbfK_a in this matrix ODE system for the acoustic pressure:\n\nmathbfM_a mathbfddotp\n+ mathbfK_a mathbfp =\nmathbff\n\nnote: Note\nThe bilinear-form function  bilform_diffusion is used to compute the matrix.\n\n\n\n\n\nacousticstiffness(self::FEMMAcoustNICE, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the acoustic mass matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler\ngeom = geometry field\nP = acoustic (perturbation) pressure field\n\nReturn a matrix.\n\n\n\n\n\n","category":"function"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustModule.acousticmass","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustModule.acousticmass","text":"acousticmass(self::FEMMAcoust, assembler::A,\n  geom::NodalField,\n  Pddot::NodalField{T}) where {T<:Number,\n  A<:AbstractSysmatAssembler}\n\nCompute the acoustic mass matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler\ngeom = geometry field\nPddot = second order rate of the acoustic (perturbation) pressure field\n\nThe acoustic \"mass\" matrix is by convention called mass, however its mechanical meaning is quite different. (It has to do with potential energy.) It is the matrix mathbfM_a in this matrix ODE system for the acoustic pressure:\n\nmathbfM_a mathbfddotp\n+ mathbfK_a mathbfp =\nmathbff\n\nnote: Note\nThe bilinear-form function  bilform_dot is used to compute the matrix.\n\n\n\n\n\n","category":"function"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustModule.inspectintegpoints","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustModule.inspectintegpoints","text":"inspectintegpoints(self::FEMMAcoust,\n    geom::NodalField{GFT},\n    P::NodalField{T},\n    temp::NodalField{FT},\n    felist::VecOrMat{IntT},\n    inspector::F,\n    idat,\n    quantity = :gradient;\n    context...) where {T <: Number, GFT, FT, IntT, F <: Function}\n\nInspect integration point quantities.\n\nArguments\n\ngeom - reference geometry field\nP - pressure field\ntemp - temperature field (ignored)\nfelist - indexes of the finite elements that are to be inspected:   The fes to be included are: fes[felist].\ncontext - struct: see the update!() method of the material.\ninspector - function with the signature       idat = inspector(idat, j, conn, x, out, loc);  where   idat - a structure or an array that the inspector may          use to maintain some state,  for instance gradient, j is the         element number, conn is the element connectivity, out is the         output of the update!() method,  loc is the location of the         integration point in the reference configuration.\n\nOutput\n\nThe updated inspector data is returned.\n\n\n\n\n\n","category":"function"},{"location":"man/man.html#Acoustics:-surface","page":"Manual","title":"Acoustics: surface","text":"","category":"section"},{"location":"man/man.html","page":"Manual","title":"Manual","text":"FEMMAcoustSurf\nacousticABC\nacousticrobin\npressure2resultantforce\npressure2resultanttorque\nacousticcouplingpanels","category":"page"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.FEMMAcoustSurf","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.FEMMAcoustSurf","text":"FEMMAcoustSurf{ID<:IntegDomain, M} <: AbstractFEMM\n\nClass for linear acoustics finite element modeling machine.\n\n\n\n\n\n","category":"type"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticABC","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticABC","text":"acousticABC(self::FEMMAcoustSurf, assembler::A,\n  geom::NodalField,\n  Pdot::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the acoustic ABC (Absorbing Boundary Condition) matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler; must be able to assemble unsymmetric matrix\ngeom = geometry field\nPdot = rate of the acoustic (perturbation) pressure field\n\nWe assume here that the impedance of this boundary is ho c.\n\n\n\n\n\n","category":"function"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticrobin","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticrobin","text":"acousticrobin(\n    self::FEMMAcoustSurf,\n    assembler::A,\n    geom::NodalField,\n    Pdot::NodalField{T},\n    impedance\n) where {T<:Number,A<:AbstractSysmatAssembler}\n\nCompute the acoustic \"Robin boundary condition\" (damping) matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler; must be able to assemble unsymmetric matrix\ngeom = geometry field\nPdot = rate of the acoustic (perturbation) pressure field\nimpedance = acoustic impedance of the boundary\n\nWe assume here that the impedance of this boundary is ho c.\n\n\n\n\n\n","category":"function"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultantforce","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultantforce","text":"pressure2resultantforce(self::FEMMAcoustSurf, assembler::A,\n  geom::NodalField,\n  P::NodalField{T},\n   Force::Field) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the rectangular coupling matrix that transcribes given pressure on the surface into the resultant force acting on the surface.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler; must be able to assemble unsymmetric matrix\ngeom = geometry field\nP = acoustic (perturbation) pressure field\nForce = field for the force resultant\n\n\n\n\n\n","category":"function"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultanttorque","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultanttorque","text":"pressure2resultanttorque(self::FEMMAcoustSurf, assembler::A,\n  geom::NodalField,\n  P::NodalField{T},\n  Torque::GeneralField, CG::FFltVec) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the rectangular coupling matrix that transcribes given pressure on the surface into the resultant torque acting on the surface with respect to the CG.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler; must be able to assemble unsymmetric matrix\ngeom = geometry field\nP = acoustic (perturbation) pressure field\nTorque = field for the torque resultant\n\n\n\n\n\n","category":"function"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticcouplingpanels","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticcouplingpanels","text":"acousticcouplingpanels(self::FEMMAcoustSurf, geom::NodalField, u::NodalField{T}) where {T}\n\nCompute the acoustic pressure-structure coupling matrix.\n\nThe acoustic pressure-nodal force matrix transforms the pressure distributed along the surface to forces acting on the nodes of the finite element model. Its transpose transforms displacements (or velocities, or accelerations) into the normal component of the displacement (or velocity, or acceleration) along the surface.\n\nArguments\n\ngeom=geometry field\nu = displacement field\n\nnote: Note\n\n\nn = outer normal (pointing into the acoustic medium).\nThe pressures along the surface are assumed constant (uniform) along   each finite element –- panel. The panel pressures are assumed to be   given the same numbers as the serial numbers of the finite elements in   the set.\n\n\n\n\n\n","category":"function"},{"location":"man/man.html#Algorithms","page":"Manual","title":"Algorithms","text":"","category":"section"},{"location":"man/man.html#Acoustics","page":"Manual","title":"Acoustics","text":"","category":"section"},{"location":"man/man.html","page":"Manual","title":"Manual","text":"AlgoAcoustModule.steadystate","category":"page"},{"location":"man/man.html#FinEtoolsAcoustics.AlgoAcoustModule.steadystate","page":"Manual","title":"FinEtoolsAcoustics.AlgoAcoustModule.steadystate","text":"steadystate(modeldata::FDataDict)\n\nSteady-state acoustics solver.\n\nmodeldata = dictionary with string keys\n\n\"fens\"  = finite element node set\n\"regions\"  = array of region dictionaries\n\"essential_bcs\" = array of essential boundary condition dictionaries\n\"ABCs\" = array of absorbing boundary condition dictionaries\n\"flux_bcs\" = array of flux boundary condition dictionaries\n\nFor each region (connected piece of the domain made of a particular material), mandatory, the  region dictionary  contains items:\n\n\"femm\" = finite element mmodel machine (mandatory);\n\nFor essential boundary conditions (optional) each dictionary would hold\n\n\"pressure\" = fixed (prescribed) pressure (scalar),  or         a function with signature             function T = f(x)         If not given, zero pressure assumed.\n\"node_list\" = list of nodes on the boundary to which the condition applies         (mandatory)\n\nFor absorbing boundary conditions (optional) each dictionary may hold\n\n\"femm\" = finite element mmodel machine (mandatory).\n\nFor flux boundary conditions (optional) each dictionary would hold\n\n\"femm\" = finite element mmodel machine (mandatory);\n\"normal_flux\" = normal component of the flux through the boundary (scalar),     which is the normal derivative of the pressure.\n\nOutput\n\nmodeldata = the dictionary is augmented with\n\n\"geom\" = the nodal field that is the geometry\n\"P\" = the nodal field that is the computed pressure (in the general a           complex-number field)\n\n\n\n\n\n","category":"function"},{"location":"man/man.html#Material-models-for-acoustics","page":"Manual","title":"Material models for acoustics","text":"","category":"section"},{"location":"man/man.html","page":"Manual","title":"Manual","text":"MatAcoustFluid\nMatAcoustFluidModule.bulkmodulus","category":"page"},{"location":"man/man.html#FinEtoolsAcoustics.MatAcoustFluidModule.MatAcoustFluid","page":"Manual","title":"FinEtoolsAcoustics.MatAcoustFluidModule.MatAcoustFluid","text":"MatAcoustFluid <: AbstractMat\n\nType for acoustic fluid material.\n\n\n\n\n\n","category":"type"},{"location":"man/man.html#FinEtoolsAcoustics.MatAcoustFluidModule.bulkmodulus","page":"Manual","title":"FinEtoolsAcoustics.MatAcoustFluidModule.bulkmodulus","text":"bulkmodulus(self::MatAcoustFluid)\n\nReturn the bulk modulus.\n\n\n\n\n\n","category":"function"},{"location":"man/man.html#Modules","page":"Manual","title":"Modules","text":"","category":"section"},{"location":"man/man.html","page":"Manual","title":"Manual","text":"FinEtoolsAcoustics.FinEtoolsAcoustics\nFinEtoolsAcoustics.FEMMAcoustModule\nFinEtoolsAcoustics.FEMMAcoustSurfModule\nFinEtoolsAcoustics.FEMMAcoustNICE\nFinEtoolsAcoustics.AlgoAcoustModule\nFinEtoolsAcoustics.MatAcoustFluidModule","category":"page"},{"location":"man/man.html#FinEtoolsAcoustics.FinEtoolsAcoustics","page":"Manual","title":"FinEtoolsAcoustics.FinEtoolsAcoustics","text":"FinEtools (C) 2017-2024, Petr Krysl\n\nFinite Element tools.  Julia implementation  of the finite element method for continuum mechanics. Package for acoustic problems.\n\n\n\n\n\n","category":"module"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustModule","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustModule","text":"FEMMAcoustModule\n\nModule for operations on interiors of domains to construct system matrices and system vectors for linear acoustics.\n\n\n\n\n\n","category":"module"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustSurfModule","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule","text":"FEMMAcoustSurfModule\n\nModule for operations on boundaries of domains to construct system matrices and system vectors for linear acoustics.\n\n\n\n\n\n","category":"module"},{"location":"man/man.html#FinEtoolsAcoustics.FEMMAcoustNICEModule.FEMMAcoustNICE","page":"Manual","title":"FinEtoolsAcoustics.FEMMAcoustNICEModule.FEMMAcoustNICE","text":"FEMMAcoustNICE{S<:AbstractFESet, F<:Function, M} <: AbstractFEMM\n\nType for linear acoustics finite element modeling machine.\n\n\n\n\n\n","category":"type"},{"location":"man/man.html#FinEtoolsAcoustics.AlgoAcoustModule","page":"Manual","title":"FinEtoolsAcoustics.AlgoAcoustModule","text":"AlgoAcoustModule\n\nModule for linear acoustics algorithms.\n\n\n\n\n\n","category":"module"},{"location":"man/man.html#FinEtoolsAcoustics.MatAcoustFluidModule","page":"Manual","title":"FinEtoolsAcoustics.MatAcoustFluidModule","text":"MatAcoustFluidModule\n\nModule for acoustic-fluid  material.\n\n\n\n\n\n","category":"module"}]
}
