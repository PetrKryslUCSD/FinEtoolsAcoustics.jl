"""
    FEMMAcoustSurfModule

Module for operations on boundaries of domains to construct system matrices and
system vectors for linear acoustics.
"""
module FEMMAcoustSurfModule

import Base.Complex

using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.FENodeSetModule: FENodeSet
using FinEtools.FESetModule: AbstractFESet, gradN!, nodesperelem, manifdim
using ..MatAcoustFluidModule: MatAcoustFluid, bulkmodulus
using FinEtools.MatModule: massdensity
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobiansurface
using FinEtools.FieldModule: ndofs, gatherdofnums!, gathervalues_asmat!, nalldofs
using FinEtools.NodalFieldModule: NodalField
using FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
using FinEtools.GeneralFieldModule: GeneralField
using FinEtools.AssemblyModule:
    AbstractSysvecAssembler,
    AbstractSysmatAssembler,
    SysmatAssemblerSparseSymm,
    startassembly!,
    assemble!,
    makematrix!,
    SysmatAssemblerSparse
using FinEtools.FEMMBaseModule: AbstractFEMM, bilform_dot
using FinEtools.MatrixUtilityModule:
    add_mggt_ut_only!, add_nnt_ut_only!, complete_lt!, locjac!
using LinearAlgebra: norm, cross
using FinEtools.DataCacheModule: DataCache

"""
    FEMMAcoustSurf{ID<:IntegDomain, M} <: AbstractFEMM

Class for linear acoustics finite element modeling machine.
"""
mutable struct FEMMAcoustSurf{ID<:IntegDomain,M} <: AbstractFEMM
    integdomain::ID # geometry data
    material::M # material object
end

"""
    acousticABC(self::FEMMAcoustSurf, assembler::A,
      geom::NodalField,
      Pdot::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}

Compute the acoustic ABC (Absorbing Boundary Condition) matrix.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler; must be able to assemble unsymmetric matrix
- `geom` = geometry field
- `Pdot` = rate of the acoustic (perturbation) pressure field

We assume here that the impedance of this boundary is ``\rho c``.
"""
function acousticABC(
    self::FEMMAcoustSurf,
    assembler::A,
    geom::NodalField,
    Pdot::NodalField{T},
) where {T<:Number,A<:AbstractSysmatAssembler}
    bulk_modulus = bulkmodulus(self.material)
    mass_density = massdensity(self.material)
    c = sqrt(bulk_modulus / mass_density) # sound speed
    oc = 1.0 / c
    return bilform_dot(self, assembler, geom, Pdot, DataCache(oc); m = 2)
end

function acousticABC(
    self::FEMMAcoustSurf,
    geom::NodalField,
    Pdot::NodalField{T},
) where {T<:Number}
    # Were we supplied assembler object?  If not make a default.
    assembler = SysmatAssemblerSparseSymm()
    return acousticABC(self, assembler, geom, Pdot)
end


"""
    acousticrobin(
        self::FEMMAcoustSurf,
        assembler::A,
        geom::NodalField,
        Pdot::NodalField{T},
        impedance
    ) where {T<:Number,A<:AbstractSysmatAssembler}

Compute the acoustic "Robin boundary condition" (damping) matrix.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler; must be able to assemble unsymmetric matrix
- `geom` = geometry field
- `Pdot` = rate of the acoustic (perturbation) pressure field
- `impedance` = acoustic impedance of the boundary

We assume here that the impedance of this boundary is ``\rho c``.
"""
function acousticrobin(
    self::FEMMAcoustSurf,
    assembler::A,
    geom::NodalField,
    Pdot::NodalField{T},
    impedance
) where {T<:Number,A<:AbstractSysmatAssembler}
    mass_density = massdensity(self.material)
    Y = mass_density / impedance
    return bilform_dot(self, assembler, geom, Pdot, DataCache(Y); m = 2)
end

function acousticrobin(
    self::FEMMAcoustSurf,
    geom::NodalField,
    Pdot::NodalField{T},
    impedance
) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm()
    return acousticrobin(self, assembler, geom, Pdot, impedance)
end

"""
    pressure2resultantforce(self::FEMMAcoustSurf, assembler::A,
      geom::NodalField,
      P::NodalField{T},
       Force::Field) where {T<:Number, A<:AbstractSysmatAssembler}

Compute the rectangular coupling matrix that transcribes given pressure
on the surface into the resultant force acting on the surface.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler; must be able to assemble unsymmetric matrix
- `geom` = geometry field
- `P` = acoustic (perturbation) pressure field
- `Force` = field for the force resultant
"""
function pressure2resultantforce(
    self::FEMMAcoustSurf,
    assembler::A,
    geom::NodalField,
    P::NodalField{T},
    Force::GeneralField,
    surfacenormal::SurfaceNormal,
) where {T<:Number,A<:AbstractSysmatAssembler}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes) # number of finite elements in the set
    ndn = ndofs(P) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    edim = ndn * nne          # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    transposedNs = [reshape(N, 1, nne) for N in Ns]
    ecoords = fill(zero(FFlt), nne, ndofs(geom)) # array of Element coordinates
    Ge = fill(zero(FFlt), 3, nne) # element coupling matrix -- used as a buffer
    coldofnums = zeros(FInt, 1, edim) # degree of freedom array -- used as a buffer
    rowdofnums = zeros(FInt, 1, 3) # degree of freedom array -- used as a buffer
    loc = fill(zero(FFlt), 1, sdim) # quadrature point location -- used as a buffer
    n = fill(zero(FFlt), 3) # normal vector -- used as a buffer
    J = fill(zero(FFlt), sdim, mdim) # Jacobian matrix -- used as a buffer
    gatherdofnums!(Force, rowdofnums, [1 2 3])# retrieve degrees of freedom
    startassembly!(assembler, 3 * edim * count(fes), 3, nalldofs(P))
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(Ge, 0.0) # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j])
            @assert Jac != 0.0
            n = updatenormal!(surfacenormal, loc, J, i, j)
            ffactor = (Jac * w[j])
            Ge = Ge + (ffactor * n) * transposedNs[j]
        end # Loop over quadrature points
        gatherdofnums!(P, coldofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, Ge, rowdofnums, coldofnums)# assemble unsymmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function pressure2resultantforce(
    self::FEMMAcoustSurf,
    geom::NodalField,
    P::NodalField{T},
    Force::GeneralField,
) where {T<:Number}
    assembler = SysmatAssemblerSparse()
    sdim = ndofs(geom)
    return pressure2resultantforce(self, assembler, geom, P, Force, SurfaceNormal(sdim))
end

function pressure2resultantforce(
    self::FEMMAcoustSurf,
    geom::NodalField,
    P::NodalField{T},
    Force::GeneralField,
    surfacenormal::SurfaceNormal,
) where {T<:Number}
    assembler = SysmatAssemblerSparse()
    return pressure2resultantforce(self, assembler, geom, P, Force, surfacenormal)
end

"""
    pressure2resultanttorque(self::FEMMAcoustSurf, assembler::A,
      geom::NodalField,
      P::NodalField{T},
      Torque::GeneralField, CG::FFltVec) where {T<:Number, A<:AbstractSysmatAssembler}

Compute the rectangular coupling matrix that transcribes given pressure
on the surface into the resultant torque acting on the surface with respect
to the CG.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler; must be able to assemble unsymmetric matrix
- `geom` = geometry field
- `P` = acoustic (perturbation) pressure field
- `Torque` = field for the torque resultant
"""
function pressure2resultanttorque(
    self::FEMMAcoustSurf,
    assembler::A,
    geom::NodalField,
    P::NodalField{T},
    Torque::GeneralField,
    CG::FFltVec,
    surfacenormal::SurfaceNormal,
) where {T<:Number,A<:AbstractSysmatAssembler}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes) # number of finite elements in the set
    ndn = ndofs(P) # number of degrees of freedom per node
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    edim = ndn * nne          # dimension of the element matrix
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    transposedNs = [reshape(N, 1, nne) for N in Ns]
    ecoords = fill(zero(FFlt), nne, ndofs(geom)) # array of Element coordinates
    Ge = fill(zero(FFlt), 3, nne) # element coupling matrix -- a buffer
    coldofnums = zeros(FInt, 1, edim) # degree of freedom array -- a buffer
    rowdofnums = zeros(FInt, 1, 3) # degree of freedom array -- a buffer
    loc = fill(zero(FFlt), 1, sdim) # quadrature point location -- a buffer
    n = fill(zero(FFlt), 3) # normal vector -- a buffer
    J = fill(zero(FFlt), sdim, mdim) # Jacobian matrix -- a buffer
    gatherdofnums!(Torque, rowdofnums, [1 2 3])# retrieve degrees of freedom
    startassembly!(assembler, 3 * edim * count(fes), 3, nalldofs(P))
    for i in eachindex(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(Ge, 0.0) # Initialize element matrix
        for j in 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j])
            n = updatenormal!(surfacenormal, loc, J, i, j)
            ffactor = (Jac * w[j])
            Ge = Ge + (ffactor * cross(vec(vec(loc) - CG), n)) * transposedNs[j]
        end # Loop over quadrature points
        gatherdofnums!(P, coldofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, Ge, rowdofnums, coldofnums)# assemble unsymmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function pressure2resultanttorque(
    self::FEMMAcoustSurf,
    geom::NodalField,
    P::NodalField{T},
    Torque::GeneralField,
    CG::FFltVec,
) where {T<:Number}
    assembler = SysmatAssemblerSparse()
    sdim = ndofs(geom)
    return pressure2resultanttorque(
        self,
        assembler,
        geom,
        P,
        Torque,
        CG,
        SurfaceNormal(sdim),
    )
end

function pressure2resultanttorque(
    self::FEMMAcoustSurf,
    geom::NodalField,
    P::NodalField{T},
    Torque::GeneralField,
    CG::FFltVec,
    surfacenormal::SurfaceNormal,
) where {T<:Number}
    assembler = SysmatAssemblerSparse()
    return pressure2resultanttorque(self, assembler, geom, P, Torque, CG, surfacenormal)
end

"""
    acousticcouplingpanels(self::FEMMAcoustSurf, geom::NodalField, u::NodalField{T}) where {T}

Compute the acoustic pressure-structure coupling matrix.

The acoustic pressure-nodal force matrix transforms the pressure distributed
along the surface to forces acting on the nodes of the finite element model.
Its transpose transforms displacements (or velocities, or accelerations) into
the normal component of the displacement (or velocity, or acceleration) along
the surface.

# Arguments
- `geom`=geometry field
- `u` = displacement field

!!! note

- `n` = outer normal (pointing into the acoustic medium).
- The pressures along the surface are assumed constant (uniform) along
    each finite element –- panel. The panel pressures are assumed to be
    given the same numbers as the serial numbers of the finite elements in
    the set.
"""
function acousticcouplingpanels(
    self::FEMMAcoustSurf,
    assembler::A,
    geom::NodalField,
    u::NodalField{T},
    surfacenormal::SurfaceNormal,
) where {A<:AbstractSysmatAssembler,T}
    fes = self.integdomain.fes
    # Constants
    nne = nodesperelem(fes) # number of nodes per element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes)     # manifold dimension of the element
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    transposedNs = [reshape(N, 1, nne) for N in Ns]
    ecoords = fill(zero(FFlt), nne, ndofs(geom)) # array of Element coordinates
    coldofnums = zeros(FInt, 1, 1) # degree of freedom array -- a buffer
    rowdofnums = zeros(FInt, 1, sdim * nne) # degree of freedom array -- a buffer
    loc = fill(zero(FFlt), 1, sdim) # quadrature point location -- a buffer
    n = fill(zero(FFlt), sdim) # normal vector -- a buffer
    J = fill(zero(FFlt), sdim, mdim) # Jacobian matrix -- a buffer
    Ge = fill(zero(FFlt), sdim * nne, 1) # Element matrix -- a buffer
    startassembly!(assembler, prod(size(Ge)) * count(fes), nalldofs(u), count(fes))
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i])
        fill!(Ge, 0.0) # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobiansurface(self.integdomain, J, loc, fes.conn[i], Ns[j])
            n = updatenormal!(surfacenormal, loc, J, i, j)
            Ge =
                Ge +
                (Jac * w[j]) *
                reshape(reshape(n, sdim, 1) * transposedNs[j], size(Ge, 1), size(Ge, 2))
        end # Loop over quadrature points
        coldofnums[1] = i
        gatherdofnums!(u, rowdofnums, fes.conn[i])# retrieve degrees of freedom
        assemble!(assembler, Ge, rowdofnums, coldofnums)# assemble unsymmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function acousticcouplingpanels(
    self::FEMMAcoustSurf,
    geom::NodalField,
    u::NodalField{T},
) where {T}
    assembler = SysmatAssemblerSparse() # The matrix is not symmetric
    sdim = ndofs(geom)
    return acousticcouplingpanels(self, assembler, geom, u, SurfaceNormal(sdim))
end

function acousticcouplingpanels(
    self::FEMMAcoustSurf,
    geom::NodalField,
    u::NodalField{T},
    surfacenormal::SurfaceNormal,
) where {T}
    assembler = SysmatAssemblerSparse() # The matrix is not symmetric
    return acousticcouplingpanels(self, assembler, geom, u, surfacenormal)
end

end
