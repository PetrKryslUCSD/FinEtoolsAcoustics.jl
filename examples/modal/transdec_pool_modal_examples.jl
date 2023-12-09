module transdec_pool_modal_examples
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection, vtkwrite
using FinEtoolsAcoustics
using DataDrop
using LinearAlgebra
using Arpack: eigs

function _run_transdec_pool(name, neigvs = 200, save_vtks = false, modelist = [])
    rho = 1000.0*phun("kg/m^3");# mass density
    c  = 1500.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    R = 0.5*phun("ft");# radius of the piston
    rim_distance = 150 * phun("ft")
    tolerance = R/100


    t0 = time()

    data = joinpath(dirname(@__FILE__()), "../data")
    input = "transdec-pool-intact-quarter.inp"
    if !isfile(joinpath(data, input))
        success(run(`unzip -qq -d $(data) $(joinpath(data, "data.zip"))`; wait = true))
    end
    mesh = import_ABAQUS(joinpath(data, input))
    fens = mesh["fens"]
    fes = mesh["fesets"][1]

    fens.xyz *= phun("ft")

    @info "$(count(fens)) nodes"

    fens1, fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 25.0, 56.65].*phun("ft"); renumb =  c ->  c[[1, 3, 2, 4]])
    fens, newfes1, fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(fens, fes, [0.0, 0.0, -1.0], [0.0, 25.0, 56.65].*phun("ft"); renumb =  c ->  c[[1, 3, 2, 4]])
    fens, newfes1, fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    bfes = meshboundary(fes)
    if save_vtks
        File =  "$(name)-boundary.vtk"
        vtkexportmesh(File, fens, bfes)
    end

    l4 = selectelem(fens, bfes, box = [-Inf, Inf, 25.0, 25.0, -Inf, Inf].*phun("ft"), inflate = tolerance)

    if save_vtks
        File =  "$(name)-free.vtk"
        vtkexportmesh(File, fens, subset(bfes, l4))
    end

    println("Pre-processing time elapsed  =  ",time() - t0,"s")

    t1  =  time()

    material = MatAcoustFluid(bulk, rho)
    # Region of the fluid
    femm = FEMMAcoust(IntegDomain(fes, TetRule(1)), material)

    geom = NodalField(fens.xyz);
    P = NodalField(zeros(nnodes(geom),1));
    setebc!(P, connectednodes(subset(bfes, l4)), true, 1, 0.0);
    applyebc!(P);
    numberdofs!(P);

    Ma = acousticmass(femm, geom, P);
    Ka = acousticstiffness(femm, geom, P);
    Ma_ff = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka_ff = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]
    d,v,nev,nconv = eigs(Ka_ff, Ma_ff; nev=neigvs, which=:SM, explicittransform=:none)
    v = real.(v)
    fs=real(sqrt.(complex(d)))./(2*pi)
    println("Eigenvalues: $fs [Hz]")

    # Visualization

    if save_vtks
        @info "Dumping visualization"
        File =  "$(name)-modes.vtk"
        scalarllist = Any[]
        for n   in modelist
            scattersysvec!(P, v[:, n])
            push!(scalarllist, ("Pressure_mode_$n", deepcopy(P.values)));
        end
        vtkexportmesh(File, fens, fes; scalars=scalarllist)
    # @async run(`"paraview.exe" $File`)
    end


    @info "Done"

    true
end #

function transdec_pool_modal(neigvs = 20)
    _run_transdec_pool("transdec_pool_modal", neigvs, false)
end

function allrun()
    println("#####################################################") 
    println("# transdec_pool_modal ")
    transdec_pool_modal()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module transdec_pool_transient_examples

nothing
