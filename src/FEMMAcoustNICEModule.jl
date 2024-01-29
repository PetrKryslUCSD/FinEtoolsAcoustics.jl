"""
    FEMMAcoustNICEModule

Acoustics with Nodally-Integrated Continuum Elements (NICE).

"""
module FEMMAcoustNICEModule

__precompile__(true)

using FinEtools.FTypesModule:
    FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.FENodeSetModule: FENodeSet
using FinEtools.FESetModule: AbstractFESet, FESetH8, FESetT4, manifdim, nodesperelem, gradN!
using FinEtools.IntegDomainModule: IntegDomain, integrationdata, Jacobianvolume
using FinEtools.FieldModule:
    ndofs,
    gatherdofnums!,
    gatherfixedvalues_asvec!,
    gathervalues_asvec!,
    gathervalues_asmat!,
    nalldofs
using FinEtools.NodalFieldModule: NodalField, nnodes
using FinEtools.FENodeToFEMapModule: FENodeToFEMap
import ..FEMMAcoustModule: acousticstiffness
using FinEtools.AssemblyModule:
    AbstractSysvecAssembler,
    AbstractSysmatAssembler,
    SysmatAssemblerSparseSymm,
    startassembly!,
    assemble!,
    makematrix!,
    makevector!,
    SysvecAssembler
using FinEtools.MatrixUtilityModule:
    add_mggt_ut_only!, complete_lt!, loc!, jac!, locjac!, adjugate3!
using FinEtools.FEMMBaseModule: AbstractFEMM
import FinEtools.FEMMBaseModule: associategeometry!
using FinEtools.MatModule: massdensity
using LinearAlgebra: mul!, Transpose, UpperTriangular, eigvals, det
using LinearAlgebra: norm, qr, diag, dot, cond, I, cross
using Statistics: mean


mutable struct _NodalBasisFunctionGradients
    gradN::FFltMat
    patchconn::FIntVec
    Vpatch::FFlt
end

"""
    FEMMAcoustNICE{S<:AbstractFESet, F<:Function, M} <: AbstractFEMM

Type for linear acoustics finite element modeling machine.
"""
mutable struct FEMMAcoustNICE{ID<:IntegDomain,M} <: AbstractFEMM
    integdomain::ID # geometry data
    material::M # material object
    associated::Bool
    nodalbasisfunctiongrad::Vector{_NodalBasisFunctionGradients}
end

function FEMMAcoustNICE(
    integdomain::IntegDomain{S,F},
    material::M,
) where {S<:FESetT4,F<:Function,M}
    return FEMMAcoustNICE(integdomain, material, false, _NodalBasisFunctionGradients[])
end

function _patchconn(fes, gl, thisnn)
    # Generate patch connectivity for a given node (thisnn)
    # from the connectivities of the finite elements attached to it.
    return vcat(
        collect(setdiff(Set([i for j = eachindex(gl) for i in fes.conn[gl[j]]]), thisnn)),
        [thisnn],
    )
end

function _buffers1(self, geom::NodalField)
    fes = self.integdomain.fes
    nne = nodesperelem(fes) # number of nodes for element
    sdim = ndofs(geom)            # number of space dimensions
    mdim = manifdim(fes) # manifold dimension of the element
    # Prepare buffers
    loc = fill(zero(FFlt), 1, sdim) # quadrature point location -- buffer
    J = fill(zero(FFlt), sdim, mdim) # Jacobian matrix -- buffer
    adjJ = fill(zero(FFlt), sdim, mdim) # Jacobian matrix -- buffer
    csmatTJ = fill(zero(FFlt), mdim, mdim) # intermediate result -- buffer
    gradN = fill(zero(FFlt), nne, mdim)
    xl = fill(zero(FFlt), nne, mdim)
    return loc, J, adjJ, csmatTJ, gradN, xl
end

function _computenodalbfungrads(self, geom)
    # # Compute the nodal basis function gradients.
    # # Return the cell array of structures with attributes
    # %      bfun_gradients{nix}.Nspd= basis function gradient matrix
    # #        bfun_gradients{nix}.Vpatch= nodal patch volume
    # #        bfun_gradients{nix}.patchconn= nodal patch connectivity

    fes = self.integdomain.fes
    npts, Ns, gradNparams, w, pc = integrationdata(self.integdomain)
    loc, J, adjJ, csmatTJ, gradN, xl = _buffers1(self, geom)

    # Get the inverse map from finite element nodes to geometric cells
    fen2fe = FENodeToFEMap(fes.conn, nnodes(geom))
    # Initialize the nodal gradients, nodal patch, and patch connectivity
    bfungrads =
        fill(_NodalBasisFunctionGradients(fill(0.0, 0, 0), fill(0, 0), 0.0), nnodes(geom))
    # Now loop over all finite element nodes in the map
    lnmap = fill(0, length(fen2fe.map)) # Local node map: buffer to speed up operations
    for nix = 1:length(fen2fe.map)
        gl = fen2fe.map[nix]
        thisnn = nix # We are at this node
        if !isempty(gl) # This node has an element patch in this block
            # establish local numbering of all nodes of the patch @ node thisnn
            p = _patchconn(fes, gl, thisnn)
            np = length(p)
            lnmap[p] .= 1:np# now store the local numbers
            c = reshape(geom.values[thisnn, :], 1, ndofs(geom))
            gradNavg = fill(0.0, np, ndofs(geom))# preallocate strain-displacement matrix
            Vpatch = 0.0
            for k in eachindex(gl)
                i = gl[k]
                kconn = collect(fes.conn[i])
                pci = findfirst(cx -> cx == thisnn, kconn)# at which node in the element are we with this quadrature point?
                @assert 1 <= pci <= nodesperelem(fes)
                # centered coordinates of nodes in the material coordinate system
                for cn in eachindex(kconn)
                    xl[cn, :] = (reshape(geom.values[kconn[cn], :], 1, ndofs(geom)) - c)
                end
                jac!(J, xl, gradNparams[pci])
                Jac = Jacobianvolume(self.integdomain, J, c, fes.conn[i], Ns[pci])
                Vpatch += Jac * w[pci]
                sgradN = gradNparams[pci] * adjugate3!(adjJ, J)
                gradNavg[lnmap[kconn], :] += (w[pci] .* sgradN)
            end
            if Vpatch == 0
                @show geom.values[thisnn, :]
                @show p
            end
            gradNavg ./= Vpatch
            bfungrads[nix] = _NodalBasisFunctionGradients(gradNavg, p, Vpatch)
            lnmap[p] .= 0 # Restore the buffer to pristine condition
        end
    end
    self.nodalbasisfunctiongrad = bfungrads
    self.associated = true
    return self
end

"""
    associategeometry!(self::F,  geom::NodalField{FFlt}) where {F<:FEMMDeforLinearESNICET4}

Associate geometry field with the FEMM.

Compute the  correction factors to account for  the shape of the  elements.
"""
function associategeometry!(self::F, geom::NodalField{FFlt}) where {F<:FEMMAcoustNICE}
    fes = self.integdomain.fes
    return _computenodalbfungrads(self, geom)
end

"""
    acousticstiffness(self::FEMMAcoustNICE, assembler::A, geom::NodalField, P::NodalField{T}) where {T<:Number, A<:AbstractSysmatAssembler}

Compute the acoustic mass matrix.

# Arguments
- `self`   =  acoustics model
- `assembler`  =  matrix assembler
- `geom` = geometry field
- `P` = acoustic (perturbation) pressure field

Return a matrix.
"""
function acousticstiffness(
    self::FEMMAcoustNICE,
    assembler::A,
    geom::NodalField,
    P::NodalField{T},
) where {T<:Number,A<:AbstractSysmatAssembler}
    @assert self.associated
    fes = self.integdomain.fes
    elmatsizeguess = 4 * nodesperelem(fes) * ndofs(P)
    startassembly!(assembler, elmatsizeguess^2 * nnodes(P), nalldofs(P), nalldofs(P))
    for nix = 1:length(self.nodalbasisfunctiongrad)
        gradN = self.nodalbasisfunctiongrad[nix].gradN
        patchconn = self.nodalbasisfunctiongrad[nix].patchconn
        Vpatch = self.nodalbasisfunctiongrad[nix].Vpatch
        nd = length(patchconn) * ndofs(P)
        elmat = fill(0.0, nd, nd) # Can we SPEED it UP?
        add_mggt_ut_only!(elmat, gradN, Vpatch)
        complete_lt!(elmat)
        dofnums = fill(0, nd)
        gatherdofnums!(P, dofnums, patchconn) # retrieve degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums) # assemble symmetric matrix
    end # Loop over elements
    return makematrix!(assembler)
end

function acousticstiffness(
    self::FEMMAcoustNICE,
    geom::NodalField,
    P::NodalField{T},
) where {T<:Number}
    # Make the default assembler object.
    assembler = SysmatAssemblerSparseSymm()
    return acousticstiffness(self, assembler, geom, P)
end

end
