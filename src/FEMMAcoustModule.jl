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
using FinEtools.FEMMBaseModule: AbstractFEMM
using FinEtools.MatrixUtilityModule: add_mggt_ut_only!, add_nnt_ut_only!, complete_lt!, locjac!
using LinearAlgebra: norm

"""
    FEMMAcoust{S<:AbstractFESet, F<:Function, M} <: AbstractFEMM

Type for linear acoustics finite element modeling machine.
"""
mutable struct FEMMAcoust{ID<:IntegDomain, M} <: AbstractFEMM
    integdomain::ID # geometry data
    material::M # material object
end

function  buffers(self::FEMMAcoust, geom::NodalField{FFlt}, P::NodalField{F}) where {F}
    fes = self.integdomain.fes
    # Constants
    nfes = count(fes); # number of finite elements in the set
    ndn = ndofs(P); # number of degrees of freedom per node
    nne =  nodesperelem(fes); # number of nodes per element
    sdim =  ndofs(geom);            # number of space dimensions
    mdim = manifdim(fes);     # manifold dimension of the element
    Cedim = ndn*nne;          # dimension of the element matrix
    # Prepare assembler and temporaries
    ecoords = fill(zero(FFlt), nne, ndofs(geom)); # array of Element coordinates
    dofnums = fill(zero(FInt), Cedim); # degree of freedom array -- used as a buffer
    loc = fill(zero(FFlt), 1, sdim); # quadrature point location -- used as a buffer
    J = fill(zero(FFlt), sdim, mdim); # Jacobian matrix -- used as a buffer
    gradN = fill(zero(FFlt), nne, mdim); # intermediate result -- used as a buffer
    elmat = fill(zero(FFlt), Cedim, Cedim);       # element matrix -- used as a buffer
    elvec = fill(zero(F), Cedim); # buffer
    elvecfix = fill(zero(F), Cedim); # buffer
    return ecoords, dofnums, loc, J, gradN, elmat, elvec, elvecfix
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
"""
function acousticmass(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}
    fes = self.integdomain.fes
    ecoords, dofnums, loc, J, gradN, elmat, elvec, elvecfix =   buffers(self, geom, P)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    Jac = 0.0;
    afactor = T(0.0);
    startassembly!(assembler, prod(size(elmat)) * count(fes),
        nalldofs(P), nalldofs(P));
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i]);
        fill!(elmat, T(0.0)); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            # gradient WRT global Cartesian coordinates
            gradN!(fes, gradN, gradNparams[j], J);
            afactor = (Jac*w[j]);
            add_mggt_ut_only!(elmat, gradN, afactor)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(P, dofnums, fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function acousticmass(self::FEMMAcoust, geom::NodalField, P::NodalField{T}) where {T<:Number}
    # Make the default assembler object.
    assembler  =  SysmatAssemblerSparseSymm();
    return acousticmass(self, assembler, geom, P);
end

"""
    nzebcloadsacousticmass(self::FEMMAcoust, assembler::A,
      geom::NodalField, P::NodalField{T}) where {T<:Number,
      A<:AbstractSysvecAssembler}

Compute load vector for nonzero EBC for prescribed pressure.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler
- `geom` = geometry field
- `P` = acoustic (perturbation) pressure field
"""
function nzebcloadsacousticmass(self::FEMMAcoust, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysvecAssembler}
    fes = self.integdomain.fes
    ecoords, dofnums, loc, J, gradN, elmat, elvec, elvecfix = buffers(self, geom, P)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    startassembly!(assembler, nalldofs(P));
    # Now loop over all finite elements in the set
    for i = 1:count(fes) # Loop over elements
        gatherfixedvalues_asvec!(P, elvecfix, fes.conn[i]);# retrieve element coordinates
        if norm(elvecfix, Inf) !=  0.0     # Is the load nonzero?
            gathervalues_asmat!(geom, ecoords, fes.conn[i]);
            fill!(elmat, T(0.0));
            for j = 1:npts # Loop over quadrature points
                locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
                Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
                # gradient WRT global Cartesian coordinates
                gradN!(fes, gradN, gradNparams[j], J);
                afactor = (Jac*w[j]);
                add_mggt_ut_only!(elmat, gradN, afactor)
            end # Loop over quadrature points
            complete_lt!(elmat)
            gatherdofnums!(P, dofnums, fes.conn[i]);# retrieve degrees of freedom
            assemble!(assembler, -elmat*elvecfix, dofnums); # assemble element load vector
        end
    end
    return  makevector!(assembler);
end

function nzebcloadsacousticmass(self::FEMMAcoust, geom::NodalField, P::NodalField{T}) where {T<:Number}
    assembler  =  SysvecAssembler(P.values[1])
    return nzebcloadsacousticmass(self, assembler, geom, P);
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
    fes = self.integdomain.fes
    ecoords, dofnums, loc, J, gradN, elmat, elvec, elvecfix = buffers(self, geom, Pddot)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    # Material
    bulk_modulus  =  bulkmodulus(self.material);
    mass_density  =  massdensity(self.material);
    c  =  sqrt(bulk_modulus/mass_density); # sound speed
    oc2 = 1.0/c^2;
    startassembly!(assembler, prod(size(elmat)) * count(fes),
        nalldofs(Pddot), nalldofs(Pddot));
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom, ecoords, fes.conn[i]);
        fill!(elmat, T(0.0)); # Initialize element matrix
        for j = 1:npts # Loop over quadrature points
            locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
            Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
            ffactor = Jac*oc2*w[j]
            add_nnt_ut_only!(elmat, Ns[j], ffactor)
        end # Loop over quadrature points
        complete_lt!(elmat)
        gatherdofnums!(Pddot, dofnums, fes.conn[i]);# retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums);# assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler);
end

function acousticstiffness(self::FEMMAcoust, geom::NodalField, Pddot::NodalField{T}) where {T<:Number}
    # Make the default assembler object.
    assembler  =  SysmatAssemblerSparseSymm();
    return acousticstiffness(self, assembler, geom, Pddot);
end

"""
    nzebcloadsacousticstiffness(self::FEMMAcoust, assembler::A,
      geom::NodalField,
      Pddot::NodalField{T}) where {T<:Number,
      A<:AbstractSysvecAssembler}

Compute load vector for nonzero EBC for prescribed second-order pressure rate.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler
- `geom` = geometry field
- `Pddot` = second order rate of the acoustic (perturbation) pressure field
"""
function nzebcloadsacousticstiffness(self::FEMMAcoust, assembler::A, geom::NodalField, Pddot::NodalField{T}) where {T<:Number, A<:AbstractSysvecAssembler}
    fes = self.integdomain.fes
    ecoords, dofnums, loc, J, gradN, elmat, elvec, elvecfix = buffers(self, geom, Pddot)
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    npts, Ns, gradNparams, w, pc  =  integrationdata(self.integdomain);
    # Material
    bulk_modulus  =   bulkmodulus(self.material);
    mass_density  =   massdensity(self.material);
    c  =  sqrt(bulk_modulus/mass_density); # sound speed
    oc2 = 1.0/c^2;
    startassembly!(assembler, nalldofs(Pddot));
    # Now loop over all finite elements in the set
    for i = 1:count(fes) # Loop over elements
        gatherfixedvalues_asvec!(Pddot, elvecfix, fes.conn[i]);# retrieve element coordinates
        if norm(elvecfix, Inf) !=  0.0  # Is the load nonzero?
            gathervalues_asmat!(geom, ecoords, fes.conn[i]);
            fill!(elmat, T(0.0));
            for j = 1:npts # Loop over quadrature points
                locjac!(loc, J, ecoords, Ns[j], gradNparams[j])
                Jac = Jacobianvolume(self.integdomain, J, loc, fes.conn[i], Ns[j]);
                ffactor = Jac*oc2*w[j]
                add_nnt_ut_only!(elmat, Ns[j], ffactor)
            end # Loop over quadrature points
            complete_lt!(elmat)
            gatherdofnums!(Pddot, dofnums, fes.conn[i]); # retrieve degrees of freedom
            assemble!(assembler, -elmat*elvecfix, dofnums); # assemble element load vector
        end
    end
    return makevector!(assembler);
end

function nzebcloadsacousticstiffness(self::FEMMAcoust, geom::NodalField, Pddot::NodalField{T}) where {T<:Number}
    assembler  =  SysvecAssembler(T(0.0))#
    return nzebcloadsacousticstiffness(self, assembler, geom, Pddot)
end

end
