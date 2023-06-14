"""
    FEMMAcoustModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear acoustics.
"""
module FEMMAcoustModule

import Base.Complex

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.FENodeSetModule: FENodeSet
using FinEtools.FESetModule: AbstractFESet, gradN!, nodesperelem, manifdim
using ..MatAcoustFluidModule: MatAcoustFluid, bulkmodulus
using FinEtools.MatModule: massdensity
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
using FinEtools.FieldModule: ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asmat!, nalldofs
using FinEtools.NodalFieldModule: NodalField
using FinEtools.AssemblyModule: AbstractSysvecAssembler, AbstractSysmatAssembler, SysmatAssemblerSparseSymm, startassembly!, assemble!, makematrix!, SysvecAssembler, makevector!
using FinEtools.FEMMBaseModule: AbstractFEMM, bilform_diffusion, bilform_dot
using FinEtools.MatrixUtilityModule: add_mggt_ut_only!, add_nnt_ut_only!, complete_lt!, locjac!
using LinearAlgebra: norm
using FinEtools.DataCacheModule: DataCache

"""
    FEMMAcoust{S<:AbstractFESet, F<:Function, M} <: AbstractFEMM

Type for linear acoustics finite element modeling machine.
"""
mutable struct FEMMAcoust{ID<:IntegDomain, M} <: AbstractFEMM
    integdomain::ID # geometry data
    material::M # material object
end

"""
    acousticmass(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}

Compute the acoustic mass matrix.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler
- `geom` = geometry field
- `P` = acoustic (perturbation) pressure field

Return a matrix.

!!! note

    The bilinear-form function  [`bilform_diffusion`](@ref) is used to compute
    the matrix.
"""
function acousticmass(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}
    f = DataCache(1.0)
    return bilform_diffusion(self, assembler, geom, P, f)
end

function acousticmass(self::FEMMAcoust, geom::NodalField, P::NodalField{T}) where {T<:Number}
    # Make the default assembler object.
    assembler  =  SysmatAssemblerSparseSymm();
    return acousticmass(self, assembler, geom, P);
end

"""
    acousticstiffness(self::FEMMAcoust, assembler::A,
      geom::NodalField,
      Pddot::NodalField{T}) where {T<:Number,
      A<:AbstractSysmatAssembler}

Compute the acoustic stiffness matrix.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler
- `geom` = geometry field
- `Pddot` = second order rate of the acoustic (perturbation) pressure field
"""
function acousticstiffness(self::FEMMAcoust, assembler::A, geom::NodalField, Pddot::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}
    # Material
    bulk_modulus  =  bulkmodulus(self.material);
    mass_density  =  massdensity(self.material);
    c  =  sqrt(bulk_modulus/mass_density); # sound speed
    oc2 = 1.0/c^2;
    f = DataCache(oc2)
    return bilform_dot(self, assembler, geom, Pddot, f)
end

function acousticstiffness(self::FEMMAcoust, geom::NodalField, Pddot::NodalField{T}) where {T<:Number}
    # Make the default assembler object.
    assembler  =  SysmatAssemblerSparseSymm();
    return acousticstiffness(self, assembler, geom, Pddot);
end

end
