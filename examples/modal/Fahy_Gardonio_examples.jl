module Fahy_Gardonio_examples
using FinEtools
using FinEtoolsAcoustics
using PlotlyLight
import Arpack: eigs

function h20_mesh()
    Lx = 0.414 * phun("m") # dimension
    Ly = 0.314 * phun("m") # dimension
    Lz = 0.360 * phun("m") # dimension
    n = 16
    fens, fes = H8block(Lx, Ly, Lz, n, n, n) # Mesh
    fens, fes = H8toH20(fens, fes)
    fens, fes
end

function h8_mesh()
    Lx = 0.414 * phun("m") # dimension
    Ly = 0.314 * phun("m") # dimension
    Lz = 0.360 * phun("m") # dimension
    n = 16
    fens, fes = H8block(Lx, Ly, Lz, n, n, n) # Mesh
    fens, fes
end

const meshfun = Dict("h20" => h20_mesh, "h8" => h8_mesh)

function fahy_gardonio_example(name)
    rho = 1.21 * phun("kg/m^3") # mass density
    c = 343.0 * phun("m/s")  # sound speed
    bulk = c^2 * rho
    neigvs = 9
    OmegaShift = 1.0

    println(
        """
        Example from Sound and Structural Vibration, Radiation, Transmission
        and Response: Second Edition, Frank J. Fahy, Paolo Gardonio, page 496.

        $(name) mesh.
        """,
        )

    t0 = time()

    fens, fes = meshfun[name]()

    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz, 1), 1))

    numberdofs!(P)

    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 3)), MatAcoustFluid(bulk, rho))

    Ma = acousticmass(femm, geom, P)
    Ka = acousticstiffness(femm, geom, P)

    d, v, nev, nconv = eigs(Ka + OmegaShift * Ma, Ma; nev = neigvs, which = :SM)
    d = d .- OmegaShift
    v = real.(v)
    fs = real(sqrt.(complex(d))) / (2 * pi)
    println("Eigenvalues: $fs [Hz]")


    println("Total time elapsed = ", time() - t0, "s")

    File = "$(name).vtk"
    en = 9
    vtkexportmesh(
        File,
        fens,
        fes,
        scalars = [("Pressure_mode_$en", v[:, en])],
    )
    # @async run(`"paraview.exe" $File`)
    println("Done")
    true

end #

function allrun()
    println("#####################################################")
    println("# fahy_gardonio_example(\"h20\") ")
    fahy_gardonio_example("h20")
    println("#####################################################")
    println("# fahy_gardonio_example(\"h8\") ")
    fahy_gardonio_example("h8")
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module Fahy_examples
nothing
