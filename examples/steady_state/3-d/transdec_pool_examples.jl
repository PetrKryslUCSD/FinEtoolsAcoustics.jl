module transdec_pool_examples
using FinEtools
using FinEtoolsAcoustics
using LinearAlgebra

function transdec_pool_example(freq = 100)
    rho = 1000.0 * phun("kg/m^3")# mass density
    c = 1500.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    omega = freq * phun("rev/s")      # frequency of the piston
    a_piston = 0.1 * omega^2 * phun("mm/s^2")     # amplitude of the piston acceleration
    R = 0.5 * phun("ft")# radius of the piston
    tolerance = R / 100


    t0 = time()

    data = joinpath(dirname(@__FILE__()), "../data")
    input = "transdec-pool-quarter.inp"
    if !isfile(joinpath(data, input))
        success(run(`unzip -qq -d $(data) $(joinpath(data, "data.zip"))`; wait = true))
    end
    mesh = import_ABAQUS(joinpath(data, input))
    fens = mesh["fens"]
    fes = mesh["fesets"][1]

    fens.xyz *= phun("ft")

    @info "$(count(fens)) nodes"

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
    File = "transdec-pool-boundary.vtk"
    vtkexportmesh(File, fens, bfes)
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

    File = "transdec-pool-pos.vtk"
    vtkexportmesh(File, fens, subset(bfes, lpos))

    File = "transdec-pool-neg.vtk"
    vtkexportmesh(File, fens, subset(bfes, lneg))

    File = "transdec-pool-free.vtk"
    vtkexportmesh(File, fens, subset(bfes, l4))

    # Piston surface mesh
    piston_pos_fes = subset(bfes, lpos)
    piston_neg_fes = subset(bfes, lneg)

    println("Pre-processing time elapsed  =  ", time() - t0, "s")

    t1 = time()

    material = MatAcoustFluid(bulk, rho)
    # Region of the fluid
    region1 = FDataDict("femm" => FEMMAcoust(IntegDomain(fes, TetRule(1)), material))

    # Surface of the piston
    flux_pos = FDataDict(
        "femm" => FEMMAcoustSurf(IntegDomain(piston_pos_fes, TriRule(3)), material),
        "normal_flux" => 1.0im * rho * a_piston,
    )
    flux_neg = FDataDict(
        "femm" => FEMMAcoustSurf(IntegDomain(piston_neg_fes, TriRule(3)), material),
        "normal_flux" => 1.0im * rho * a_piston,
    )

    ebc1 =
        FDataDict("node_list" => connectednodes(subset(bfes, l4)), "pressure" => x -> 0.0)

    # Make model data
    modeldata = FDataDict(
        "fens" => fens,
        "omega" => omega,
        "regions" => [region1],
        "essential_bcs" => [ebc1],
        "flux_bcs" => [flux_pos, flux_neg],
    )

    # Call the solver
    modeldata = FinEtoolsAcoustics.AlgoAcoustModule.steadystate(modeldata)

    println("Computing time elapsed  =  ", time() - t1, "s")
    println("Total time elapsed  =  ", time() - t0, "s")

    geom = modeldata["geom"]
    P = modeldata["P"]

    File = "transdec-pool-f=$(freq)-p.vtk"
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
    @async run(`"paraview.exe" $File`)
end #

function allrun()
    println("#####################################################")
    println("# transdec_pool_example ")
    transdec_pool_example(100)
    transdec_pool_example(200)
    transdec_pool_example(500)
    transdec_pool_example(750)
    transdec_pool_example(1000)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module baffled_piston_examples
nothing
