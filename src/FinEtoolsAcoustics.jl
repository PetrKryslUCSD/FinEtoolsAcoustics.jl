"""
FinEtools (C) 2017-2019, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics. Package for heat diffusion problems.
"""
module FinEtoolsAcoustics

__precompile__(true)

include("allmodules.jl")

# Exports follow:

###########################################################################
# Acoustics functionality
###########################################################################
using FinEtoolsAcoustics.MatAcoustFluidModule: MatAcoustFluid
# Exported: type of acoustic fluid material
export MatAcoustFluid

using FinEtoolsAcoustics.FEMMAcoustModule: FEMMAcoust, acousticmass, nzebcloadsacousticmass, acousticstiffness, nzebcloadsacousticstiffness
# Exported: type for linear acoustics  and discretization methods
export FEMMAcoust, acousticmass, nzebcloadsacousticmass,
    acousticstiffness, nzebcloadsacousticstiffness

using FinEtoolsAcoustics.FEMMAcoustNICEModule: FEMMAcoustNICE, acousticmass
# Exported: type for linear acoustics  and discretization methods
export FEMMAcoustNICE, acousticmass

using FinEtoolsAcoustics.FEMMAcoustSurfModule: FEMMAcoustSurf, acousticABC, pressure2resultantforce, pressure2resultanttorque, acousticcouplingpanels
# Exported: type for acoustic absorbing boundary condition  and  transformation matrices from pressure  to resultants
export FEMMAcoustSurf, acousticABC, pressure2resultantforce, pressure2resultanttorque, acousticcouplingpanels

end # module
