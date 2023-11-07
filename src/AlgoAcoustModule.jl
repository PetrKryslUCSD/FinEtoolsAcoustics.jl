"""
    AlgoAcoustModule

Module for linear acoustics algorithms.
"""
module AlgoAcoustModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.AlgoBaseModule: dcheck!, matrix_blocked, vector_blocked
using FinEtools.FieldModule: ndofs, setebc!, numberdofs!, applyebc!, scattersysvec!, nfreedofs, nalldofs, gathersysvec
using FinEtools.NodalFieldModule: NodalField, nnodes
using FinEtools.FEMMBaseModule: associategeometry!, distribloads
using ..FEMMAcoustModule: acousticmass, acousticstiffness
using ..FEMMAcoustSurfModule: acousticABC
using FinEtools.ForceIntensityModule: ForceIntensity
using SparseArrays: spzeros
using LinearAlgebra: norm, lu, cross

"""
    steadystate(modeldata::FDataDict)

Steady-state acoustics solver.

`modeldata` = dictionary with string keys

- `"fens"`  = finite element node set
- `"regions"`  = array of region dictionaries
- `"essential_bcs"` = array of essential boundary condition dictionaries
- `"ABCs"` = array of absorbing boundary condition dictionaries
- `"flux_bcs"` = array of flux boundary condition dictionaries

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains items:
- `"femm"` = finite element mmodel machine (mandatory);

For essential boundary conditions (optional) each dictionary
would hold
  + `"pressure"` = fixed (prescribed) pressure (scalar),  or
            a function with signature
                function T = f(x)
            If not given, zero pressure assumed.
  + `"node_list"` = list of nodes on the boundary to which the condition applies
            (mandatory)

For absorbing boundary conditions (optional) each dictionary
may hold
  + `"femm"` = finite element mmodel machine (mandatory).

For flux boundary conditions (optional) each dictionary
would hold
  + `"femm"` = finite element mmodel machine (mandatory);
  + `"normal_flux"` = normal component of the flux through the boundary (scalar),
        which is the normal derivative of the pressure.

# Output
`modeldata` = the dictionary is augmented with
- `"geom"` = the nodal field that is the geometry
- `"P"` = the nodal field that is the computed pressure (in the general a
            complex-number field)
"""
function steadystate(modeldata::FDataDict)
    # Lists of recognized keys for the data dictionaries:
    modeldata_recognized_keys = ["fens", "regions", "essential_bcs", "flux_bcs"]
    essential_bcs_recognized_keys = ["pressure", "node_list"]
    ABCs_recognized_keys = ["femm"]
    flux_bcs_recognized_keys = ["femm", "normal_flux"]
    regions_recognized_keys = ["femm"]

    # Check that keys defined  in the dictionary are matched
    dcheck!(modeldata, modeldata_recognized_keys)

    # Extract the nodes
    fens = get(()->error("Must get fens!"), modeldata, "fens")

    # Construct the geometry field
    geom = NodalField(fens.xyz)

    # Construct the acoustic pressure field
    P = NodalField(zeros(FCplxFlt,size(fens.xyz,1),1))

    # The following properties are the same for all regions and all boundaries (ABC)
    omega = get(()->error("Must get angular frequency!"), modeldata, "omega");

    # Apply the essential boundary conditions on the acoustic pressure field
    essential_bcs = get(modeldata, "essential_bcs", nothing);
    if (essential_bcs !== nothing)
        for j in eachindex(essential_bcs)
            ess = essential_bcs[j]
            dcheck!(ess, essential_bcs_recognized_keys)
            fenids = get(()->error("Must get node list!"), ess, "node_list");
            pressure = get(ess, "pressure", nothing);
            Pfixed = zeros(FCplxFlt,length(fenids)); # default is zero pressure
            if (pressure !== nothing)
                if (typeof(pressure) <: Function) # pressure supplied through function
                    for k = 1:length(fenids)
                        Pfixed[k] = pressure(geom.values[fenids[k],:])[1];
                    end
                else # pressure given as  a constant
                    fill!(Pfixed, pressure);
                end
            end
            setebc!(P, fenids[:], true, 1, Pfixed);
            applyebc!(P);
        end
    end

    # Number the equations
    numberdofs!(P)           #,Renumbering_options); # NOT DONE <<<<<<<<<<<<<<<<<

    # Initialize the acoustic load vector
    F = gathersysvec(P, :a)

    # Construct the system acoustic mass and stiffness matrix
    # and the absorbing boundary condition (ABC) matrix
    C = spzeros(nalldofs(P), nalldofs(P)); # (all zeros, for the moment)
    S = spzeros(nalldofs(P), nalldofs(P)); # (all zeros, for the moment)
    D = spzeros(nalldofs(P), nalldofs(P)); # (all zeros, for the moment)
    regions = get(()->error("Must get regions!"), modeldata, "regions")
    for i in eachindex(regions)
        region = regions[i]
        dcheck!(region, regions_recognized_keys)
        femm = get(()->error("Must get femm for the region!"), region, "femm");
        # Add up all the acoustic mass matrices for all the regions
        C = C + acousticmass(femm, geom, P);
        S  = S + acousticstiffness(femm, geom, P);
    end

    # Compute the ABC matrices
    ABCs = get(modeldata, "ABCs", nothing);
    if (ABCs != nothing)
        for j = 1:length(ABCs)
            ABC = ABCs[j]
            dcheck!(ABC, ABCs_recognized_keys)
            femm = get(()->error("Must get femm for the ABC!"), ABC, "femm");
            D = D + acousticABC(femm, geom, P);
        end
    end

    # Process the flux boundary condition:
    # dP/dn=-rho*a ... The normal derivative of the pressure in terms of the acceleration
    flux_bcs = get(modeldata, "flux_bcs", nothing);
    if (flux_bcs != nothing)
        for j = 1:length(flux_bcs)
            fluxbc = flux_bcs[j]
            dcheck!(fluxbc, flux_bcs_recognized_keys)
            normal_flux = get(()->error("Must get normal flux value!"), fluxbc, "normal_flux");
            if (typeof(normal_flux) <: Function)
                fi = ForceIntensity(FCplxFlt, 1, normal_flux);
            else
                if typeof(normal_flux) <: AbstractArray
                else
                    normal_flux = FCplxFlt[normal_flux]
                end
                fi = ForceIntensity(normal_flux);
            end
            femm = get(()->error("Must get femm for the flux BC!"), fluxbc, "femm");
            F = F + distribloads(femm, geom, P, fi, 2);
        end
    end

    # Separate the blocks of the system matrices
    S_ff, S_fd = matrix_blocked(S, nfreedofs(P), nfreedofs(P))[(:ff, :fd)]
    C_ff, C_fd = matrix_blocked(C, nfreedofs(P), nfreedofs(P))[(:ff, :fd)]
    D_ff, D_fd = matrix_blocked(D, nfreedofs(P), nfreedofs(P))[(:ff, :fd)]
    F_f, F_d = vector_blocked(F, nfreedofs(P))[(:f, :d)]

    # Loads due to the essential boundary conditions on the pressure field
    essential_bcs = get(modeldata, "essential_bcs", nothing);
    if (essential_bcs != nothing)
        F_f = F_f - C_fd * F_d - S_fd * F_d
    end

    # Solve for the pressures
    K = lu((-omega^2*S_ff +omega*1.0im*D_ff + C_ff));
    vP = K\F_f;
    scattersysvec!(P, vP[:])

    # Update the model data
    setindex!(modeldata, geom, "geom");
    setindex!(modeldata, P, "P");
    return modeldata            # ... And return the updated model data
end

end
