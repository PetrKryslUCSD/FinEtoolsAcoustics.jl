"""
    FEMMAcoustModule

Module for operations on interiors of domains to construct system matrices and
system vectors for linear acoustics.
"""
module FEMMAcoustModule

import Base.Complex

using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.FENodeSetModule: FENodeSet
using FinEtools.FESetModule: AbstractFESet, gradN!, nodesperelem, manifdim
using ..MatAcoustFluidModule: MatAcoustFluid, bulkmodulus
using FinEtools.MatModule: massdensity
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
using FinEtools.FieldModule:
    ndofs, gatherdofnums!, gatherfixedvalues_asvec!, gathervalues_asmat!, nalldofs, gathervalues_asvec!
using FinEtools.NodalFieldModule: NodalField
using FinEtools.AssemblyModule:
    AbstractSysvecAssembler,
    AbstractSysmatAssembler,
    SysmatAssemblerSparseSymm,
    startassembly!,
    assemble!,
    makematrix!,
    SysvecAssembler,
    makevector!
using FinEtools.FEMMBaseModule: AbstractFEMM, bilform_diffusion, bilform_dot, finite_elements
using FinEtools.MatrixUtilityModule:
    add_mggt_ut_only!, add_nnt_ut_only!, complete_lt!, locjac!
using LinearAlgebra: norm, mul!
using FinEtools.DataCacheModule: DataCache

"""
    FEMMAcoust{ID<:IntegDomain, M} <: AbstractFEMM

Type for linear acoustics finite element modeling machine.
"""
mutable struct FEMMAcoust{ID<:IntegDomain,M} <: AbstractFEMM
    integdomain::ID # geometry data
    material::M # material object
end

"""
    acousticstiffness(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}

Compute the acoustic "stiffness" matrix.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler
- `geom` = geometry field
- `P` = acoustic (perturbation) pressure field

Return a matrix.

The acoustic "stiffness" matrix is by convention called stiffness, however its
mechanical meaning is quite different. (It has to do with kinetic energy.) It
is the matrix ``\\mathbf{K}_a`` in this matrix ODE system for the acoustic
pressure:

```math
\\mathbf{M}_a \\mathbf{\\ddot{p}}
+ \\mathbf{K}_a \\mathbf{{p}} =
\\mathbf{{f}}
```

!!! note

    The bilinear-form function  [`bilform_diffusion`](@ref) is used to compute
    the matrix.
"""
function acousticstiffness(
    self::FEMMAcoust,
    assembler::A,
    geom::NodalField,
    P::NodalField{T},
) where {T<:Number,A<:AbstractSysmatAssembler}
    f = DataCache(1.0)
    return bilform_diffusion(self, assembler, geom, P, f)
end

function acousticstiffness(
    self::FEMMAcoust,
    geom::NodalField,
    P::NodalField{T},
) where {T<:Number}
    # Make the default assembler object.
    assembler = SysmatAssemblerSparseSymm()
    return acousticstiffness(self, assembler, geom, P)
end

"""
    acousticmass(self::FEMMAcoust, assembler::A,
      geom::NodalField,
      Pddot::NodalField{T}) where {T<:Number,
      A<:AbstractSysmatAssembler}

Compute the acoustic mass matrix.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler
- `geom` = geometry field
- `Pddot` = second order rate of the acoustic (perturbation) pressure field

The acoustic "mass" matrix is by convention called mass, however its
mechanical meaning is quite different. (It has to do with potential energy.) It
is the matrix ``\\mathbf{M}_a`` in this matrix ODE system for the acoustic
pressure:

```math
\\mathbf{M}_a \\mathbf{\\ddot{p}}
+ \\mathbf{K}_a \\mathbf{{p}} =
\\mathbf{{f}}
```

!!! note

    The bilinear-form function  [`bilform_dot`](@ref) is used to compute
    the matrix.
"""
function acousticmass(
    self::FEMMAcoust,
    assembler::A,
    geom::NodalField,
    Pddot::NodalField{T},
) where {T<:Number,A<:AbstractSysmatAssembler}
    # Material
    bulk_modulus = bulkmodulus(self.material)
    mass_density = massdensity(self.material)
    c = sqrt(bulk_modulus / mass_density) # sound speed
    oc2 = 1.0 / c^2
    f = DataCache(oc2)
    return bilform_dot(self, assembler, geom, Pddot, f)
end

function acousticmass(
    self::FEMMAcoust,
    geom::NodalField,
    Pddot::NodalField{T},
) where {T<:Number}
    # Make the default assembler object.
    assembler = SysmatAssemblerSparseSymm()
    return acousticmass(self, assembler, geom, Pddot)
end

function _buffers1(self::FEMMAcoust,
    geom::NodalField{GFT},
    P::NodalField{FT}) where {GFT, FT}
    # Constants
    fes = finite_elements(self)
    IntT = eltype(P.dofnums)
    nfes = count(fes) # number of finite elements in the set
    ndn = ndofs(P) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)   # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    Kedim = ndn * nne      # dimension of the element matrix
    ecoords = fill(zero(GFT), nne, ndofs(geom)) # array of Element coordinates
    dofnums = fill(zero(IntT), Kedim) # buffer
    loc = fill(zero(GFT), 1, sdim) # buffer
    J = fill(zero(GFT), sdim, mdim) # buffer
    gradN = fill(zero(GFT), nne, mdim) # buffer
    Pe = fill(zero(FT), nodesperelem(fes)) # nodal pressures -- buffer
    qpgradP = fill(zero(FT), 1, sdim) # Pressure gradient -- buffer
    return ecoords,
    dofnums,
    loc,
    J,
    gradN,
    Pe,
    qpgradP
end

"""
    inspectintegpoints(self::FEMMAcoust,
        geom::NodalField{GFT},
        P::NodalField{T},
        temp::NodalField{FT},
        felist::VecOrMat{IntT},
        inspector::F,
        idat,
        quantity = :gradient;
        context...) where {T <: Number, GFT, FT, IntT, F <: Function}

Inspect integration point quantities.

# Arguments
- `geom` - reference geometry field
- `P` - pressure field
- `temp` - temperature field (ignored)
- `felist` - indexes of the finite elements that are to be inspected:
    The fes to be included are: `fes[felist]`.
- `context` - struct: see the update!() method of the material.
- `inspector` - function with the signature
        `idat = inspector(idat, j, conn, x, out, loc);`
   where
    `idat` - a structure or an array that the inspector may
           use to maintain some state,  for instance gradient, `j` is the
          element number, `conn` is the element connectivity, `out` is the
          output of the `update!()` method,  `loc` is the location of the
          integration point in the *reference* configuration.
# Output
The updated inspector data is returned.
"""
function inspectintegpoints(self::FEMMAcoust,
    geom::NodalField{GFT},
    P::NodalField{T},
    temp::NodalField{FT},
    felist::VecOrMat{IntT},
    inspector::F,
    idat,
    quantity = :gradient;
    context...) where {T <: Number, GFT, FT, IntT, F <: Function}
    fes = finite_elements(self)
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    ecoords, dofnums, loc, J, gradN, Pe, qpgradP = _buffers1(self, geom, P)
    t = 0.0
    dt = 0.0
    # Loop over  all the elements and all the quadrature points within them
    for i in felist # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        gathervalues_asvec!(P, Pe, fes.conn[i])# retrieve element temperatures
        for j in 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j])
            gradN!(fes, gradN, gradNparams[j], J)
            # Quadrature point quantities
            mul!(qpgradP, reshape(Pe, 1, :), gradN) # temperature gradient in material coordinates
            # Call the inspector
            idat = inspector(idat, i, fes.conn[i], ecoords, qpgradP, loc)
        end # Loop over quadrature points
    end # Loop over elements
    return idat # return the updated inspector data
end


end
