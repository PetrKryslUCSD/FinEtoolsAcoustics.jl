"""
Sphere modes.

For the data 
rho = 1.2*phun("kg/m^3");# mass density
c  = 340.0*phun("m/s");# sound speed
bulk =  c^2*rho;
R = 1000.0*phun("mm");# radius of the piston

the reference 

@article{GAO2013914,
title = {Eigenvalue analysis for acoustic problem in 3D by boundary element method with the block Sakurai–Sugiura method},
journal = {Engineering Analysis with Boundary Elements},
volume = {37},
number = {6},
pages = {914-923},
year = {2013},
issn = {0955-7997},
doi = {https://doi.org/10.1016/j.enganabound.2013.03.015},
url = {https://www.sciencedirect.com/science/article/pii/S0955799713000714},
author = {Haifeng Gao and Toshiro Matsumoto and Toru Takahashi and Hiroshi Isakari},
keywords = {Eigenvalues, Acoustic, The block SS method, Boundary element method, Burton–Miller's method},
abstract = {This paper presents accurate numerical solutions for nonlinear eigenvalue analysis of three-dimensional acoustic cavities by boundary element method (BEM). To solve the nonlinear eigenvalue problem (NEP) formulated by BEM, we employ a contour integral method, called block Sakurai–Sugiura (SS) method, by which the NEP is converted to a standard linear eigenvalue problem and the dimension of eigenspace is reduced. The block version adopted in present work can also extract eigenvalues whose multiplicity is larger than one, but for the complex connected region which includes a internal closed boundary, the methodology yields fictitious eigenvalues. The application of the technique is demonstrated through the eigenvalue calculation of sphere with unique homogenous boundary conditions, cube with mixed boundary conditions and a complex connected region formed by cubic boundary and spherical boundary, however, the fictitious eigenvalues can be identified by Burton–Miller's method. These numerical results are supported by appropriate convergence study and comparisons with close form.}
}

shows the wave numbers in Table 1:

The multiplicity of the Dirichlet eigenvalues.
Wavenumber*R                    Multiplicity
3.14159, 6.28319, 9.42478       1
4.49340, 7.72525, 10.90412      3
5.76346, 9.09501, 12.32294      5
6.98793, 10.41711, 13.69802     7
8.18256, 11.70491, 15.03966     9
9.35581, 12.96653, 16.35471    11

"""

module sphere_mode_examples
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtoolsAcoustics
using FinEtools.MeshExportModule
using LinearAlgebra
using Arpack: eigs

function sphere_h8_in_air()
    rho = 1.02 * phun("kg/m^3")# mass density
    c = 340.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    R = 1000.0 * phun("mm")# radius of the sphere

    nPerradius = 10# number of elements along the radius of the scatterer
    tolerance = R / nPerradius / 100
    neigvs = 20

    println("""
    
    Sphere with Dirichlet boundary conditions: modal analysis.
    Hexahedral H8 mesh.
    Exact fundamental frequency: $(c/2/R)

    The multiplicity of the Dirichlet eigenvalues.
    Wavenumber*R                    Multiplicity
    3.14159, 6.28319, 9.42478       1
    4.49340, 7.72525, 10.90412      3
    5.76346, 9.09501, 12.32294      5
    6.98793, 10.41711, 13.69802     7
    8.18256, 11.70491, 15.03966     9
    9.35581, 12.96653, 16.35471    11
    """)

    wn_table = [
        ([3.14159, 6.28319, 9.42478], 1),
        ([4.49340, 7.72525, 10.90412], 3),
        ([5.76346, 9.09501, 12.32294], 5),
        ([6.98793, 10.41711, 13.69802], 7),
        ([8.18256, 11.70491, 15.03966], 9),
        ([9.35581, 12.96653, 16.35471], 11),
    ]
    @info "Reference frequencies"
    for i in axes(wn_table, 1)
        fq = wn_table[i][1] ./ R .* c / (2 * pi)
        @info "$(fq), multiplicity $(wn_table[i][2])"
    end

    # Hexahedral mesh
    fens, fes = H8spheren(R, nPerradius)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0];
        renumb = c -> c[[1, 4, 3, 2, 5, 8, 7, 6]],
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0];
        renumb = c -> c[[1, 4, 3, 2, 5, 8, 7, 6]],
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 0.0];
        renumb = c -> c[[1, 4, 3, 2, 5, 8, 7, 6]],
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz, 1), 1))
    bfes = meshboundary(fes)
    setebc!(P, connectednodes(bfes))
    numberdofs!(P)

    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 3)), MatAcoustFluid(bulk, rho))

    Ma = acousticmass(femm, geom, P)
    Ka = acousticstiffness(femm, geom, P)
    Ma_ff = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka_ff = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]
    d, v, nconv = eigs(Ka_ff, Ma_ff; nev = neigvs, which = :SM, explicittransform = :none)

    v = real.(v)
    fs = real(sqrt.(complex(d))) ./ (2 * pi)
    @info("Frequencies: $fs [Hz]")
    @info("Wavenumbers: $((2*pi).*fs./c./phun("""m""")) [Hz]")

    File = "sphere_h8_in_air.vtk"
    scalarllist = Any[]
    for n = 1:15
        scattersysvec!(P, v[:, n])
        push!(scalarllist, ("Pressure_mode_$n", deepcopy(P.values)))
    end
    vtkexportmesh(
        File,
        connasarray(fes),
        geom.values,
        FinEtools.MeshExportModule.VTK.H8;
        scalars = scalarllist,
    )
    @async run(`"paraview.exe" $File`)

    true

end # sphere_h8_in_air

function sphere_t4_in_water()

    rho = 1000 * phun("kg/m^3")# mass density
    c = 1500.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    R = 500.0 * phun("mm")# radius of the sphere

    neigvs = 70

    println("""
    
    Sphere with Dirichlet boundary conditions: model analysis.
    Sphere of radius $(R), in WATER.
    Tetrahedral T4 mesh.
    Exact fundamental frequency: $(c/2/R)

    """)

    wn_table = [
        ([3.14159, 6.28319, 9.42478], 1),
        ([4.49340, 7.72525, 10.90412], 3),
        ([5.76346, 9.09501, 12.32294], 5),
        ([6.98793, 10.41711, 13.69802], 7),
        ([8.18256, 11.70491, 15.03966], 9),
        ([9.35581, 12.96653, 16.35471], 11),
    ]
    @info "Reference frequencies"
    for i in axes(wn_table, 1)
        fq = wn_table[i][1] ./ R .* c / (2 * pi)
        @info "$(fq), multiplicity $(wn_table[i][2])"
    end

    data = joinpath(dirname(@__FILE__()), "../data")
    input = "sphere-0_03.inp"
    if !isfile(joinpath(data, input))
        success(run(`unzip -qq -d $(data) $(joinpath(data, "data.zip"))`; wait = true))
    end
    mesh = import_ABAQUS(joinpath(data, input))
    fens = mesh["fens"]
    fes = mesh["fesets"][1]

    @info "$(count(fens)) nodes"

    t0 = time()

    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz, 1), 1))
    bfes = meshboundary(fes)
    setebc!(P, connectednodes(bfes))
    numberdofs!(P)

    femm = FEMMAcoust(IntegDomain(fes, TetRule(1)), MatAcoustFluid(bulk, rho))

    Ma = acousticmass(femm, geom, P)
    Ka = acousticstiffness(femm, geom, P)
    Ma_ff = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka_ff = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]
    d, v, nev, nconv =
        eigs(Ka_ff, Ma_ff; nev = neigvs, which = :SM, explicittransform = :none)
    v = real.(v)
    fs = real(sqrt.(complex(d))) ./ (2 * pi)
    @info("Frequencies: $fs [Hz]")
    ks = (2 * pi) .* fs ./ c ./ phun("m")
    @info("Wavenumbers: $(ks) [m]")

    @info "Time: $(time() - t0)"

    File = "sphere_t4_in_water.vtk"
    scalarllist = Any[]
    for n = 1:66
        scattersysvec!(P, v[:, n])
        push!(scalarllist, ("Pressure_mode_$n", deepcopy(P.values)))
    end
    vtkexportmesh(
        File,
        connasarray(fes),
        geom.values,
        FinEtools.MeshExportModule.VTK.T4;
        scalars = scalarllist,
    )
    @async run(`"paraview.exe" $File`)

    true

end # sphere_h8_in_air

function allrun()
    println("#####################################################")
    println("# sphere_h8_in_air ")
    sphere_h8_in_air()
    println("#####################################################")
    println("# sphere_t4_in_water ")
    @time sphere_t4_in_water()
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module sphere_mode_examples
nothing
