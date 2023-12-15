module transdec_pool_examples
using FinEtools
using FinEtoolsAcoustics
using LinearAlgebra
using DataDrop
using Plots

const R = 0.5 * phun("ft")# radius of the piston
const tolerance = R / 100

function mesh()
    center = [0.0, 6.0, 56.65] .* phun("ft")
    cog = center + [1.0, 0.0, 0.0] .* phun("m")
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
        from = center,
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

    skull_subset = selectelem(
        fens,
        fes,
        distance = 0.4 * phun("m"),
        from = cog,
        allin = false,
        inflate = tolerance,
    )

    # File = "transdec-pool-pos.vtk"
    # vtkexportmesh(File, fens, subset(bfes, lpos))

    # File = "transdec-pool-neg.vtk"
    # vtkexportmesh(File, fens, subset(bfes, lneg))

    # File = "transdec-pool-free.vtk"
    # vtkexportmesh(File, fens, subset(bfes, l4))

    File = "transdec-pool-cog.vtk"
    vtkexportmesh(File, fens, subset(fes, skull_subset))

    # Piston surface mesh
    piston_pos_fes = subset(bfes, lpos)
    piston_neg_fes = subset(bfes, lneg)

    free_surface_fes = subset(bfes, l4)

    return fens, fes, piston_pos_fes, piston_neg_fes, free_surface_fes, skull_subset
end

function transdec_pool_example(freq = [100])
    rho = 1002.0 * phun("kg/m^3")# mass density
    c = 1490.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    v_piston = 1000 * phun("Pa")  / (rho * c)   # amplitude of the piston acceleration
    # Definition of the transect
    from = [-0.01524000021336, 1.8288000000000002, 17.2669204572]
    to = from .+ [-7.5, 0, 0]
    npoints = 50
    geometricaltolerance = tolerance / 100

    # println("Pre-processing time elapsed  =  ", time() - t0, "s")

    fens, fes, piston_pos_fes, piston_neg_fes, free_surface_fes, skull_subset = mesh()

    material = MatAcoustFluid(bulk, rho)
    # Region of the fluid
    region1 = FDataDict("femm" => FEMMAcoust(IntegDomain(fes, TetRule(1)), material))


    ebc1 = FDataDict("node_list" => connectednodes(free_surface_fes), "pressure" => x -> 0.0)

    qpgrads = zeros(ComplexF64, count(fes), 3)
    for f in freq
        @info "Frequency $f"
        omega = 2 * pi * f
        # Surface of the piston
        flux_pos = FDataDict(
            "femm" => FEMMAcoustSurf(IntegDomain(piston_pos_fes, TriRule(3)), material),
            "normal_flux" => 1.0im * rho * v_piston * omega,
            )
        flux_neg = FDataDict(
            "femm" => FEMMAcoustSurf(IntegDomain(piston_neg_fes, TriRule(3)), material),
            "normal_flux" => 1.0im * rho * v_piston * omega,
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

        qpgrads .= zero(eltype(P.values))
        idat = inspectintegpoints(region1["femm"], geom, P,
            NodalField([1.0]), collect(1:count(fes)),
            (idat, i, conn, xe, out, loc) -> let
                qpgrads = idat
                qpgrads[i, :] .= out[:]
                return qpgrads
            end, qpgrads, :gradient)
        qpgrads = idat ./ region1["femm"].integdomain.integration_rule.npts

        # The integral is only over the neighborhood of the skull; everything else is blanked out
        volf = ElementalField(zeros(size(qpgrads, 1), 1))
        volf.values[skull_subset, :] .= 1.0
        gradf = ElementalField(zeros(eltype(P.values), size(qpgrads)))
        gradf.values[skull_subset, :] .= qpgrads[skull_subset, :]

        vol = integratefieldfunction(region1["femm"], geom, volf,
            (x, a) -> a[1]; initial = zero(eltype(gradf.values)), m = 3)
        g1 = integratefieldfunction(region1["femm"], geom, gradf,
            (x, a) -> a[1]; initial = zero(eltype(gradf.values)), m = 3)
        g2 = integratefieldfunction(region1["femm"], geom, gradf,
            (x, a) -> a[2]; initial = zero(eltype(gradf.values)), m = 3)
        g3 = integratefieldfunction(region1["femm"], geom, gradf,
            (x, a) -> a[3]; initial = zero(eltype(gradf.values)), m = 3)
        File = "transdec-pool-f=$(f)-pgrad.json"
        DataDrop.store_json(File, Dict("grad" => [g1, g2, g3] ./ vol))

        # File = "transdec-pool-f=$(f)-p.vtk"
        # vtkexportmesh(
        #     File,
        #     fens,
        #     fes;
        #     scalars = [
        #     ("absP", abs.(P.values)),
        #     ("realP", real.(P.values)),
        #     ("imagP", imag.(P.values)),
        #     ],
        #     )
    # @async run(`"paraview.exe" $File`)
    end
    true
end #

function allrun()
    println("#####################################################")
    println("# transdec_pool_example ")
    # transdec_pool_example(100)
    # transdec_pool_example(200)
    # transdec_pool_example(500)
    # transdec_pool_example(750)
    # transdec_pool_example(1000)
    freq = collect(170:2:1000)
    transdec_pool_example(freq)

    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module
nothing
