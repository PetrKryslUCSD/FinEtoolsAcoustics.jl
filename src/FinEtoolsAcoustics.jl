"""
FinEtools (C) 2017-2024, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics. Package for acoustic problems.
"""
module FinEtoolsAcoustics

__precompile__(true)

include("allmodules.jl")

# Enable LSP look-up in test modules.
if false include("../test/runtests.jl") end


using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
export FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict

# Exports follow:

###########################################################################
# Acoustics functionality
###########################################################################
using FinEtoolsAcoustics.MatAcoustFluidModule: MatAcoustFluid
# Exported: type of acoustic fluid material
export MatAcoustFluid

using FinEtoolsAcoustics.FEMMAcoustModule: FEMMAcoust, acousticstiffness, acousticmass, inspectintegpoints
# Exported: type for linear acoustics  and discretization methods
export FEMMAcoust, acousticstiffness, acousticmass, inspectintegpoints

using FinEtoolsAcoustics.FEMMAcoustNICEModule: FEMMAcoustNICE, acousticstiffness
# Exported: type for linear acoustics  and discretization methods
export FEMMAcoustNICE, acousticstiffness

using FinEtoolsAcoustics.FEMMAcoustSurfModule:
    FEMMAcoustSurf,
    acousticABC,
    acousticrobin,
    pressure2resultantforce,
    pressure2resultanttorque,
    acousticcouplingpanels
# Exported: type for acoustic absorbing boundary condition  and  transformation matrices from pressure  to resultants
export FEMMAcoustSurf,
    acousticABC, acousticrobin,
    pressure2resultantforce, pressure2resultanttorque, acousticcouplingpanels

end # module
