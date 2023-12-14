module transdec_pool_examples
using FinEtools
using FinEtoolsAcoustics
using LinearAlgebra
using DataDrop
using Plots

function _transect_interpolation(fens, fes, from, to, npoints, geometricaltolerance)
    parametrictolerance = 1.0e-5
    distances = Float64[]
    nodebox = initbox!([], vec(from))
    for i = 1:npoints
        p = (i - npoints) / (1 - npoints) .* from + (i - 1) / (npoints - 1) .* to
        push!(distances, norm(p - from))
        updatebox!(nodebox, vec(p))
    end
    nodebox = inflatebox!(nodebox, geometricaltolerance)
    el = selectelem(fens, fes; overlappingbox = nodebox)
    interp = []
    for i = 1:npoints
        p = (i - npoints) / (1 - npoints) .* from + (i - 1) / (npoints - 1) .* to
        for e in el
            c = [k for k in fes.conn[e]]
            pc, success = map2parametric(
                fes,
                fens.xyz[c, :],
                vec(p);
                tolerance = parametrictolerance,
                maxiter = 7,
            )
            @assert success # this shouldn't be tripped; normally we succeed
            if inparametric(fes, pc; tolerance = parametrictolerance)
                push!(interp, e)
                break
            end
        end
    end
    @assert length(interp) == npoints
    return distances, interp
end

const R = 0.5 * phun("ft")# radius of the piston
const tolerance = R / 100

function mesh()

    data = joinpath(dirname(@__FILE__()), "../../data")
    input = "transdec-pool-quarter.inp"
    if !isfile(joinpath(data, input))
        success(run(`unzip -qq -d $(data) $(joinpath(data, "data.zip"))`; wait = true))
    end
    mesh = import_ABAQUS(joinpath(data, input))
    fens = mesh["fens"]
    fes = mesh["fesets"][1]

    fens.xyz *= phun("ft")

    # @info "$(count(fens)) nodes"

    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 25.0, 56.65] .* phun("ft");
        renumb = c -> c[[1, 3, 2, 4]],
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    bfes = meshboundary(fes)
    # File = "transdec-pool-boundary.vtk"
    # vtkexportmesh(File, fens, bfes)
    # @async run(`"paraview.exe" $File`)

    l1 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])
    l2 = selectelem(fens, bfes, facing = true, direction = [+1.0 0.0 0.0])
    l3 = selectelem(
        fens,
        bfes,
        distance = R,
        from = [0.0, 6.0, 56.65] .* phun("ft"),
        inflate = tolerance,
    )
    l4 = selectelem(
        fens,
        bfes,
        box = [-Inf, Inf, 25.0, 25.0, -Inf, Inf] .* phun("ft"),
        inflate = tolerance,
    )

    lpos = intersect(l2, l3)
    lneg = intersect(l1, l3)

    # File = "transdec-pool-pos.vtk"
    # vtkexportmesh(File, fens, subset(bfes, lpos))

    # File = "transdec-pool-neg.vtk"
    # vtkexportmesh(File, fens, subset(bfes, lneg))

    # File = "transdec-pool-free.vtk"
    # vtkexportmesh(File, fens, subset(bfes, l4))

    # Piston surface mesh
    piston_pos_fes = subset(bfes, lpos)
    piston_neg_fes = subset(bfes, lneg)

    free_surface_fes = subset(bfes, l4)

    return fens, fes, piston_pos_fes, piston_neg_fes, free_surface_fes
end

function transdec_pool_example(freq = [100])
    rho = 1002.0 * phun("kg/m^3")# mass density
    c = 1490.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    d_piston = 0.1 * phun("mm")     # amplitude of the piston acceleration
    # Definition of the transect
    from = [-0.01524000021336, 1.8288000000000002, 17.2669204572]
    to = from .+ [-7.5, 0, 0]
    npoints = 50
    geometricaltolerance = tolerance / 100

    # println("Pre-processing time elapsed  =  ", time() - t0, "s")

    fens, fes, piston_pos_fes, piston_neg_fes, free_surface_fes = mesh()

    material = MatAcoustFluid(bulk, rho)
    # Region of the fluid
    region1 = FDataDict("femm" => FEMMAcoust(IntegDomain(fes, TetRule(1)), material))


    ebc1 = FDataDict("node_list" => connectednodes(free_surface_fes), "pressure" => x -> 0.0)

    function inspector(idat, i, conn, xe, out, loc)
        qpgradnorms = idat
        qpgradnorms[i] = norm(out)
        return qpgradnorms
    end

    for f in freq
        @info "Frequency $f"
        omega = 2 * pi * f
        # Surface of the piston
        flux_pos = FDataDict(
            "femm" => FEMMAcoustSurf(IntegDomain(piston_pos_fes, TriRule(3)), material),
            "normal_flux" => 1.0im * rho * d_piston * omega^2,
            )
        flux_neg = FDataDict(
            "femm" => FEMMAcoustSurf(IntegDomain(piston_neg_fes, TriRule(3)), material),
            "normal_flux" => 1.0im * rho * d_piston * omega^2,
            )
        modeldata = FDataDict(
            "fens" => fens,
            "omega" => omega,
            "regions" => [region1],
            "essential_bcs" => [ebc1],
            "flux_bcs" => [flux_pos, flux_neg],
            )
        modeldata = FinEtoolsAcoustics.AlgoAcoustModule.steadystate(modeldata)

        geom = modeldata["geom"]
        P = modeldata["P"]

        # Collect the norms of the gradients for each element
        qpgradnorms = zeros(count(fes))

        idat = inspectintegpoints(region1["femm"], geom, P, NodalField([1.0]), collect(1:count(fes)), inspector, qpgradnorms, :gradient)
        qpgradnorms = idat

        distances, interp = _transect_interpolation(fens, fes, from, to, npoints, geometricaltolerance)
        pgtrans = zeros(length(distances))
        for (j, ip) in enumerate(interp)
            pgtrans[j] = qpgradnorms[ip]
        end

        d = Dict(
            "freq" => f,
            "distances" => distances,
            "pgtrans" => pgtrans,
            )
        File = "transdec-pool-f=$(f)-pgtrans.json"
        DataDrop.store_json(File, d)

        File = "transdec-pool-f=$(f)-p.vtk"
        vtkexportmesh(
            File,
            fens,
            fes;
            scalars = [
            ("absP", abs.(P.values)),
            ("realP", real.(P.values)),
            ("imagP", imag.(P.values)),
            ],
            )
    # @async run(`"paraview.exe" $File`)
    end
    true
end #

#=
plot_transects = Main.transdec_pool_examples.plot_transects
pl = nothing
for f in 440:1:520;
   plot_transects(pl, "transdec-pool-f=$f-pgtrans.json", f);
   sleep(1)
end
    =#
function plot_transects(pl, f, freq)
    d = DataDrop.retrieve_json(f)
    pgtrans = d["pgtrans"]
    distances = d["distances"]
    gr()
    if pl == nothing
        pl = plot(
            vec(distances),
            vec(pgtrans),
            leg = false,
            ylims = (0, 1000000),
            title = "Frequency: $(freq)",
            xaxis = ("Distance along transect",),
            yaxis = ("||Pressure gradient||",),
            )
    else
        pl = plot!(pl,
            vec(distances),
            vec(pgtrans),
            leg = false,
            ylims = (0, 100000),
            title = "Frequency: $(freq)",
            xaxis = ("Distance along transect",),
            yaxis = ("||Pressure gradient||",),
            )
    end
    display(pl)
    return pl
end

function allrun()
    println("#####################################################")
    println("# transdec_pool_example ")
    # transdec_pool_example(100)
    # transdec_pool_example(200)
    # transdec_pool_example(500)
    # transdec_pool_example(750)
    # transdec_pool_example(1000)
    freq = collect(430:1:450)
    transdec_pool_example(freq)

    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module baffled_piston_examples
nothing
