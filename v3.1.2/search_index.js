var documenterSearchIndex = {"docs":
[{"location":"guide/guide.html#Guide","page":"How to guide","title":"Guide","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The FinEtools package is used here to solve heat conduction problems.","category":"page"},{"location":"guide/guide.html#Modules","page":"How to guide","title":"Modules","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The package FinEtoolsAcoustics has the following structure:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"Top-level:     FinEtoolsAcoustics is the  top-level module.  \nAcoustics: AlgoAcoustModule (algorithms), FEMMAcoustModule, FEMMAcoustSurfModule (FEM machines to evaluate the matrix and vector quantities),  MatAcoustFluidModule (acoustic fluid material).","category":"page"},{"location":"guide/guide.html#Acoustics-FEM-machines","page":"How to guide","title":"Acoustics FEM machines","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"There is one for  the interior integrals  and one for  boundary integrals. The  machine for the interior integrals can be used to compute:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"Evaluate the acoustic-mass matrix and the acoustic-stiffness matrix.\nEvaluate the load vector corresponding to prescribed pressure  or the prescribed second order  rate of the pressure.","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The machine for the boundary integrals can be used to compute:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"Compute  transformation matrix to convert  pressure  to resultant force  or pressure to resultant torque.\nCompute the acoustic  ABC  (absorbing boundary condition) matrix.","category":"page"},{"location":"guide/guide.html#Acoustics-algorithms","page":"How to guide","title":"Acoustics algorithms","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"At the moment there is one algorithm, for steady-state (harmonic) acoustics.","category":"page"},{"location":"guide/guide.html#Example:-baffled-piston","page":"How to guide","title":"Example:  baffled piston","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"After the mesh  has been generated, the modeldata can be set up: Here we begin with  the region.","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"material = MatAcoustFluid(bulk, rho)\nregion1 =  FDataDict(\"femm\"=>FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material))","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"We set up a definition of the absorbing boundary condition:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"abc1  =  FDataDict(\"femm\"=>FEMMAcoustSurf(IntegDomain(outer_fes, GaussRule(2, 2)),\n          material))","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The  surface of the piston is associated with a known-flux  boundary condition (prescribed normal component of the velocity):","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"flux1  =  FDataDict(\"femm\"=>FEMMAcoustSurf(IntegDomain(piston_fes, GaussRule(2, 2)),\n          material),  \"normal_flux\"=> -1.0im*rho*omega*v_piston);","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"And finally we make the model data,","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"modeldata =  FDataDict(\"fens\"=>  fens,\n                 \"omega\"=>omega,\n                 \"regions\"=>[region1],\n                 \"flux_bcs\"=>[flux1], \"ABCs\"=>[abc1])","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"and call  the solver:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"modeldata = FinEtools.AlgoAcoustModule.steadystate(modeldata)","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"When  the algorithm completes, modeldata[\"P\"] is the computed pressure field.","category":"page"},{"location":"index.html#FinEtoolsAcoustics-Documentation","page":"Home","title":"FinEtoolsAcoustics Documentation","text":"","category":"section"},{"location":"index.html#Tutorials","page":"Home","title":"Tutorials","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The tutorials are provided in the form of Julia scripts and Markdown files in a separate package. ","category":"page"},{"location":"index.html#Manual","page":"Home","title":"Manual","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The description of the types and the functions, organized by module and/or other logical principle.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n    \"man/types.md\",\n    \"man/functions.md\",\n]\nDepth = 2","category":"page"},{"location":"index.html#Conceptual-guide","page":"Home","title":"Conceptual guide","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The construction of the toolkit is described: the composition of modules, the basic data structures, the methodology of computing quantities required in the finite element methodology, and more.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n    \"guide/guide.md\",\n]\nDepth = 1","category":"page"},{"location":"man/types.html#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"man/types.html#FEM-machines","page":"Types","title":"FEM machines","text":"","category":"section"},{"location":"man/types.html#Acoustics","page":"Types","title":"Acoustics","text":"","category":"section"},{"location":"man/types.html","page":"Types","title":"Types","text":"Modules = [FinEtools, FinEtoolsAcoustics.FEMMAcoustModule, FinEtoolsAcoustics.FEMMAcoustSurfModule]\nPrivate = true\nOrder = [:type]","category":"page"},{"location":"man/types.html#FinEtoolsAcoustics.FEMMAcoustModule.FEMMAcoust","page":"Types","title":"FinEtoolsAcoustics.FEMMAcoustModule.FEMMAcoust","text":"FEMMAcoust{ID<:IntegDomain, M} <: AbstractFEMM\n\nType for linear acoustics finite element modeling machine.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.FEMMAcoustSurf","page":"Types","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.FEMMAcoustSurf","text":"FEMMAcoustSurf{ID<:IntegDomain, M} <: AbstractFEMM\n\nClass for linear acoustics finite element modeling machine.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#Material-models","page":"Types","title":"Material models","text":"","category":"section"},{"location":"man/types.html#Material-models-for-acoustics","page":"Types","title":"Material models for acoustics","text":"","category":"section"},{"location":"man/types.html","page":"Types","title":"Types","text":"Modules = [FinEtools, FinEtools.MatModule, FinEtoolsAcoustics.MatAcoustFluidModule]\nPrivate = true\nOrder = [:type]","category":"page"},{"location":"man/types.html#FinEtools.MatModule.AbstractMat","page":"Types","title":"FinEtools.MatModule.AbstractMat","text":"AbstractMat\n\nAbstract type of material.\n\n\n\n\n\n","category":"type"},{"location":"man/types.html#FinEtoolsAcoustics.MatAcoustFluidModule.MatAcoustFluid","page":"Types","title":"FinEtoolsAcoustics.MatAcoustFluidModule.MatAcoustFluid","text":"MatAcoustFluid <: AbstractMat\n\nType for acoustic fluid material.\n\n\n\n\n\n","category":"type"},{"location":"man/functions.html#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"man/functions.html#FEM-machines","page":"Functions","title":"FEM machines","text":"","category":"section"},{"location":"man/functions.html#Acoustics","page":"Functions","title":"Acoustics","text":"","category":"section"},{"location":"man/functions.html","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsAcoustics.FEMMAcoustModule, FinEtoolsAcoustics.FEMMAcoustSurfModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustModule.acousticmass-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoust, A, NodalField, NodalField{T}}} where {T<:Number, A<:AbstractSysmatAssembler}","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustModule.acousticmass","text":"acousticmass(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the acoustic mass matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler\ngeom = geometry field\nP = acoustic (perturbation) pressure field\n\nReturn a matrix.\n\nnote: Note\nThe bilinear-form function  bilform_diffusion is used to compute the matrix.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustModule.acousticstiffness-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoust, A, NodalField, NodalField{T}}} where {T<:Number, A<:AbstractSysmatAssembler}","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustModule.acousticstiffness","text":"acousticstiffness(self::FEMMAcoust, assembler::A,\n  geom::NodalField,\n  Pddot::NodalField{T}) where {T<:Number,\n  A<:AbstractSysmatAssembler}\n\nCompute the acoustic stiffness matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler\ngeom = geometry field\nPddot = second order rate of the acoustic (perturbation) pressure field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticABC-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoustSurf, A, NodalField, NodalField{T}}} where {T<:Number, A<:AbstractSysmatAssembler}","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticABC","text":"acousticABC(self::FEMMAcoustSurf, assembler::A,\n  geom::NodalField,\n  Pdot::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the acoustic ABC (Absorbing Boundary Condition) matrix.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler; must be able to assemble unsymmetric matrix\ngeom = geometry field\nPdot = rate of the acoustic (perturbation) pressure field\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticcouplingpanels-Union{Tuple{T}, Tuple{A}, Tuple{FEMMAcoustSurf, A, NodalField, NodalField{T}, SurfaceNormal}} where {A<:AbstractSysmatAssembler, T}","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.acousticcouplingpanels","text":"acousticcouplingpanels(self::FEMMAcoustSurf, geom::NodalField, u::NodalField{T}) where {T}\n\nCompute the acoustic pressure-structure coupling matrix.\n\nThe acoustic pressure-nodal force matrix transforms the pressure distributed along the surface to forces acting on the nodes of the finite element model. Its transpose transforms displacements (or velocities, or accelerations) into the normal component of the displacement (or velocity, or acceleration) along the surface.\n\nArguments\n\ngeom=geometry field\nu = displacement field\n\nnote: Note\n\n\nn = outer normal (pointing into the acoustic medium).\nThe pressures along the surface are assumed constant (uniform) along   each finite element –- panel. The panel pressures are assumed to be   given the same numbers as the serial numbers of the finite elements in   the set.\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultantforce-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoustSurf, A, NodalField, NodalField{T}, GeneralField, SurfaceNormal}} where {T<:Number, A<:AbstractSysmatAssembler}","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultantforce","text":"pressure2resultantforce(self::FEMMAcoustSurf, assembler::A,\n  geom::NodalField,\n  P::NodalField{T},\n   Force::Field) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the rectangular coupling matrix that transcribes given pressure on the surface into the resultant force acting on the surface.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler; must be able to assemble unsymmetric matrix\ngeom = geometry field\nP = acoustic (perturbation) pressure field\nForce = field for the force resultant\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultanttorque-Union{Tuple{A}, Tuple{T}, Tuple{FEMMAcoustSurf, A, NodalField, NodalField{T}, GeneralField, Vector{Float64}, SurfaceNormal}} where {T<:Number, A<:AbstractSysmatAssembler}","page":"Functions","title":"FinEtoolsAcoustics.FEMMAcoustSurfModule.pressure2resultanttorque","text":"pressure2resultanttorque(self::FEMMAcoustSurf, assembler::A,\n  geom::NodalField,\n  P::NodalField{T},\n  Torque::GeneralField, CG::FFltVec) where {T<:Number, A<:AbstractSysmatAssembler}\n\nCompute the rectangular coupling matrix that transcribes given pressure on the surface into the resultant torque acting on the surface with respect to the CG.\n\nArguments\n\nself   =  acoustics model\nassembler  =  matrix assembler; must be able to assemble unsymmetric matrix\ngeom = geometry field\nP = acoustic (perturbation) pressure field\nTorque = field for the torque resultant\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Algorithms","page":"Functions","title":"Algorithms","text":"","category":"section"},{"location":"man/functions.html#Acoustics-2","page":"Functions","title":"Acoustics","text":"","category":"section"},{"location":"man/functions.html","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsAcoustics.AlgoAcoustModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsAcoustics.AlgoAcoustModule.steadystate-Tuple{Dict{String, Any}}","page":"Functions","title":"FinEtoolsAcoustics.AlgoAcoustModule.steadystate","text":"steadystate(modeldata::FDataDict)\n\nSteady-state acoustics solver.\n\nmodeldata = dictionary with string keys\n\n\"fens\"  = finite element node set\n\"regions\"  = array of region dictionaries\n\"essential_bcs\" = array of essential boundary condition dictionaries\n\"ABCs\" = array of absorbing boundary condition dictionaries\n\"flux_bcs\" = array of flux boundary condition dictionaries\n\nFor each region (connected piece of the domain made of a particular material), mandatory, the  region dictionary  contains items:\n\n\"femm\" = finite element mmodel machine (mandatory);\n\nFor essential boundary conditions (optional) each dictionary would hold\n\n\"pressure\" = fixed (prescribed) pressure (scalar),  or         a function with signature             function T = f(x)         If not given, zero pressure assumed.\n\"node_list\" = list of nodes on the boundary to which the condition applies         (mandatory)\n\nFor absorbing boundary conditions (optional) each dictionary may hold\n\n\"femm\" = finite element mmodel machine (mandatory).\n\nFor flux boundary conditions (optional) each dictionary would hold\n\n\"femm\" = finite element mmodel machine (mandatory);\n\"normal_flux\" = normal component of the flux through the boundary (scalar),     which is the normal derivative of the pressure.\n\nOutput\n\nmodeldata = the dictionary is augmented with\n\n\"geom\" = the nodal field that is the geometry\n\"P\" = the nodal field that is the computed pressure (in the general a           complex-number field)\n\n\n\n\n\n","category":"method"},{"location":"man/functions.html#Material-models","page":"Functions","title":"Material models","text":"","category":"section"},{"location":"man/functions.html#Material-models-for-acoustics","page":"Functions","title":"Material models for acoustics","text":"","category":"section"},{"location":"man/functions.html","page":"Functions","title":"Functions","text":"Modules = [FinEtools, FinEtoolsAcoustics.MatAcoustFluidModule]\nPrivate = true\nOrder = [:function]","category":"page"},{"location":"man/functions.html#FinEtoolsAcoustics.MatAcoustFluidModule.bulkmodulus-Tuple{MatAcoustFluid}","page":"Functions","title":"FinEtoolsAcoustics.MatAcoustFluidModule.bulkmodulus","text":"bulkmodulus(self::MatAcoustFluid)\n\nReturn the bulk modulus.\n\n\n\n\n\n","category":"method"}]
}
