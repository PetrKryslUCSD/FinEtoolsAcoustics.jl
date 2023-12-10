module sphsld_acous_coupl_examples
using FinEtools
using FinEtoolsAcoustics
using FinEtoolsDeforLinear
using FinEtools.MeshExportModule
using Arpack: eigs
using SparseArrays
using PlotlyLight

function sphsld_acous_coupl_examples_H8()
    R = 0.5 * phun("m")
    E = 205000 * phun("MPa")
    nu = 0.3
    rho = 7850 * phun("kg/m^3")# mass density
    dummybulk, dummyrho = (1.0, 1.0)
    neigvs = 6 + 18
    OmegaShift = (2 * pi * 100.0)^2
    nperradius = 2
    tolerance = R / 10000 # geometrical tolerance

    # Construct the mesh
    r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]]
    origin = [0.0, 0.0, 0.0]
    fens, fes = H8spheren(R, nperradius) # Mesh
    fens1, fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], origin, renumb = r)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(fens, fes, [0.0, -1.0, 0.0], origin, renumb = r)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(fens, fes, [0.0, 0.0, -1.0], origin, renumb = r)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    # Extract the mesh of the boundary
    bfes = meshboundary(fes)

    # Debugging graphics
    # File  =   "Sphere.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    numberdofs!(u)

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    M = mass(femm, geom, u)

    d, v, nev, nconv = eigs(K + OmegaShift * M, M; nev = neigvs, which = :SM)
    d = d .- OmegaShift
    fs = real(sqrt.(complex(d))) / (2 * pi)
    println("Eigenvalues: $fs [Hz]")

    femm = FEMMAcoustSurf(
        IntegDomain(bfes, TrapezoidalRule(2)),
        MatAcoustFluid(dummybulk, dummyrho),
    )
    G = acousticcouplingpanels(femm, geom, u)
    I, J, V = findnz(G)
    p = PlotlyLight.Plot()
    p(x = J, y = I, mode = "markers")
    # p.layout.title.text = "Matrix G"
    p.layout.yaxis.title = "Row"
    p.layout.yaxis.range = [size(G, 1) + 1, 0]
    p.layout.xaxis.title = "Column"
    p.layout.xaxis.range = [0, size(G, 2) + 1]
    p.layout.xaxis.side = "top"
    p.layout.margin.pad = 10
    display(p)

    true

end # sphsld_acous_coupl_examples_H8

function allrun()
    println("#####################################################")
    println("# sphsld_acous_coupl_examples_H8 ")
    sphsld_acous_coupl_examples_H8()
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module sphsld_acous_coupl_examples
nothing
