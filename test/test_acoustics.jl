
module mmrrigid
using FinEtools

using FinEtoolsAcoustics
using Test
import Arpack: eigs
import LinearAlgebra: norm
function test()

    # println("""
    #
    # Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 123.
    # Internal resonance problem. Reference frequencies: 90.7895, 181.579, 215.625, 233.959, 272.368, 281.895
    # Quadrilateral mesh.
    # """)

    t0 = time()

    rho = 1.21 * 1e-9# mass density
    c = 345.0 * 1000# millimeters per second
    bulk = c^2 * rho
    Lx = 1900.0# length of the box, millimeters
    Ly = 800.0 # length of the box, millimeters
    n = 14#
    neigvs = 18
    OmegaShift = 10.0

    fens, fes = Q4block(Lx, Ly, n, n) # Mesh

    geom = NodalField(fens.xyz)
    P = NodalField(fill(zero(FFlt), size(fens.xyz, 1), 1))

    numberdofs!(P)

    femm = FEMMAcoust(IntegDomain(fes, GaussRule(2, 2)), MatAcoustFluid(bulk, rho))

    Ma = acousticmass(femm, geom, P)
    Ka = acousticstiffness(femm, geom, P)

    d, v, nev, nconv = eigs(Ka + OmegaShift * Ma, Ma; nev = neigvs, which = :SM)
    d = d .- OmegaShift
    fs = real(sqrt.(complex(d))) / (2 * pi)
    # println("Eigenvalues: $fs [Hz]")
    #
    #
    # println("Total time elapsed = ",time() - t0,"s")
    # println("$fs[2]-90.7895")
    @test abs(fs[2] - 90.9801) < 1e-4

    File = "rigid_box.vtk"
    scalarllist = Any[]
    for n = 1:15
        push!(scalarllist, ("Pressure_mode_$n", v[:, n]))
    end
    vtkexportmesh(
        File,
        connasarray(fes),
        geom.values,
        FinEtools.MeshExportModule.VTK.Q4;
        scalars = scalarllist,
    )
    sleep(1.0)
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)

    true
end
end
using .mmrrigid
mmrrigid.test()

module fahyL2example

using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsAcoustics
using Test
import Arpack: eigs
function test()
    # println("""
    # Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
    # Frank J. Fahy, Paolo Gardonio, page 483.
    #
    # 1D mesh.
    # """)

    t0 = time()

    rho = 1.21 * 1e-9# mass density
    c = 343.0 * 1000# millimeters per second
    bulk = c^2 * rho
    L = 500.0# length of the box, millimeters
    A = 200.0 # cross-sectional area of the box
    graphics = true# plot the solution as it is computed?
    n = 40#
    neigvs = 8
    OmegaShift = 10.0

    fens, fes = L2block(L, n) # Mesh

    geom = NodalField(fens.xyz)
    P = NodalField(fill(zero(FFlt), size(fens.xyz, 1), 1))

    numberdofs!(P)

    femm = FEMMAcoust(IntegDomain(fes, GaussRule(1, 2)), MatAcoustFluid(bulk, rho))

    Ma = acousticmass(femm, geom, P)
    Ma = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka = acousticstiffness(femm, geom, P)
    Ka = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]

    d, v, nev, nconv = eigs(Ka + OmegaShift * Ma, Ma; nev = neigvs, which = :SM)
    d = d .- OmegaShift
    fs = real(sqrt.(complex(d))) / (2 * pi)
    # println("Eigenvalues: $fs [Hz]")


    # println("Total time elapsed  =  ",time() - t0,"s")

    # using Plots
    # plotly()
    # en = 2
    # ix = sortperm(geom.values[:])
    # plot(geom.values[:][ix], v[:,en][ix], color = "blue",
    # title = "Fahy example, mode $en" , xlabel = "x", ylabel = "P")
    # gui()

    @test (fs[2] - 343.08) / 343.08 < 0.001

    #println("Total time elapsed = ",time() - t0,"s")

    # using Winston
    # en=2
    # pl = FramedPlot(title="Fahy example, mode $en",xlabel="x",ylabel="P")
    # setattr(pl.frame, draw_grid=true)
    # ix=sortperm(geom.values[:])
    # add(pl, Curve(geom.values[:][ix],v[:,en][ix], color="blue"))

    # display(pl)

    true
end
end
using .fahyL2example
fahyL2example.test()

module mmfahyH8example
using FinEtools

using FinEtoolsAcoustics
using Test
import Arpack: eigs
function test()

    # println("""
    # Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
    # Frank J. Fahy, Paolo Gardonio, page 483.
    #
    # Hexahedral mesh.
    # """)

    t0 = time()

    rho = 1.21 * 1e-9# mass density
    c = 343.0 * 1000# millimeters per second
    bulk = c^2 * rho
    L = 500.0# length of the box, millimeters
    A = 200.0 # cross-sectional area of the box
    n = 40#
    neigvs = 8
    OmegaShift = 10.0

    fens, fes = H8block(L, sqrt(A), sqrt(A), n, 1, 1) # Mesh

    geom = NodalField(fens.xyz)
    P = NodalField(fill(zero(FFlt), size(fens.xyz, 1), 1))

    numberdofs!(P)

    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), MatAcoustFluid(bulk, rho))


    Ma = acousticmass(femm, geom, P)
    Ka = acousticstiffness(femm, geom, P)

    d, v, nev, nconv = eigs(Ka + OmegaShift * Ma, Ma; nev = neigvs, which = :SM)
    d = d .- OmegaShift
    fs = real(sqrt.(complex(d))) / (2 * pi)
    # println("Eigenvalues: $fs [Hz]")
    #
    #
    # println("Total time elapsed = ",time() - t0,"s")

    # File =  "fahy_H8.vtk"
    # en = 5;
    # vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
    #           scalars=[("Pressure_mode_$en", v[:,en])])
    # println("Done")

    @test (fs[2] - 343.08) / 343.08 < 0.001

    true
end
end
using .mmfahyH8example
mmfahyH8example.test()

module mmfahyH27example
using FinEtools
using FinEtoolsAcoustics
using Test
import Arpack: eigs
function test()

    # println("""
    # Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
    # Frank J. Fahy, Paolo Gardonio, page 483.
    #
    # Hexahedral mesh.
    # """)

    t0 = time()

    rho = 1.21 * 1e-9# mass density
    c = 343.0 * 1000# millimeters per second
    bulk = c^2 * rho
    L = 500.0# length of the box, millimeters
    A = 200.0 # cross-sectional area of the box
    n = 40#
    neigvs = 8
    OmegaShift = 10.0

    fens, fes = H8block(L, sqrt(A), sqrt(A), n, 1, 1) # Mesh
    fens, fes = H8toH27(fens, fes)

    geom = NodalField(fens.xyz)
    P = NodalField(fill(zero(FFlt), size(fens.xyz, 1), 1))

    numberdofs!(P)

    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 3)), MatAcoustFluid(bulk, rho))


    Ma = acousticmass(femm, geom, P)
    Ka = acousticstiffness(femm, geom, P)

    d, v, nev, nconv = eigs(Ka + OmegaShift * Ma, Ma; nev = neigvs, which = :SM)
    d = d .- OmegaShift
    fs = real(sqrt.(complex(d))) / (2 * pi)
    # println("Eigenvalues: $fs [Hz]")
    #
    #
    # println("Total time elapsed = ",time() - t0,"s")

    # File =  "fahy_H8.vtk"
    # en = 5;
    # vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
    #           scalars=[("Pressure_mode_$en", v[:,en])])
    # println("Done")

    @test (fs[2] - 343.08) / 343.08 < 0.001

    true
end
end
using .mmfahyH27example
mmfahyH27example.test()

module mstraight_duct_H8_1
using FinEtools

using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsAcoustics
using Test
function test()
    t0 = time()

    rho = 1.21 * phun("kg/m^3")# mass density
    c = 343.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    omega = 54.5901 * phun("rev/s")
    vn0 = -1.0 * phun("m/s")
    Lx = 10.0 * phun("m")# length of the box, millimeters
    Ly = 1.0 * phun("m") # length of the box, millimeters
    n = 20#number of elements along the length

    # println("""
    #
    # Straight duct with anechoic termination.
    # Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 44.
    # Both real and imaginary components of the pressure should have amplitude of
    # rho*c = $(rho*c).
    #
    # Hexahedral mesh.
    # """)

    fens, fes = H8block(Lx, Ly, Ly, n, 2, 2) # Mesh
    bfes = meshboundary(fes)
    L0 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])
    L10 = selectelem(fens, bfes, facing = true, direction = [+1.0 0.0 0.0])
    nLx = selectnode(fens, box = [0.0 Lx 0.0 0.0 0.0 0.0], inflate = Lx / 1.0e5)

    geom = NodalField(fens.xyz)
    P = NodalField(fill(zero(Complex{Float64}), size(fens.xyz, 1), 1))

    numberdofs!(P)


    material = MatAcoustFluid(bulk, rho)
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material)

    Ma = acousticmass(femm, geom, P)
    Ma = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka = acousticstiffness(femm, geom, P)
    Ka = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]

    E10femm = FEMMAcoustSurf(IntegDomain(subset(bfes, L10), GaussRule(2, 2)), material)
    D = acousticABC(E10femm, geom, P)

    E0femm = FEMMBase(IntegDomain(subset(bfes, L0), GaussRule(2, 2)))
    fi = ForceIntensity(-1.0im * omega * rho * vn0)
    F = distribloads(E0femm, geom, P, fi, 2)
    F = vector_blocked(F, nfreedofs(P))

    p = (-omega^2 * Ma + omega * 1.0im * D + Ka) \ F.f
    scattersysvec!(P, p[:])

    # println("Pressure amplitude bounds: ")
    # println("  real $(minimum(real(P.values)))/$(maximum(real(P.values)))")
    # println("  imag $(minimum(imag(P.values)))/$(maximum(imag(P.values)))")
    #
    # println("Total time elapsed  =  ",time() - t0,"s")

    File = "straight_duct.vtk"
    scalars = real(P.values)
    vtkexportmesh(
        File,
        connasarray(fes),
        geom.values,
        FinEtools.MeshExportModule.VTK.H8;
        scalars = [("Pressure", scalars)],
    )
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)

    ref = rho * c
    @test abs(minimum(real(P.values)) - (-ref)) / ref < 0.02
    @test abs(minimum(imag(P.values)) - (-ref)) / ref < 0.02
    @test abs(maximum(real(P.values)) - (+ref)) / ref < 0.02
    @test abs(maximum(imag(P.values)) - (+ref)) / ref < 0.02

    # plotly()
    # ix = sortperm(geom.values[nLx,1])
    # plot(geom.values[nLx,1][ix], real(P.values)[nLx][ix], color = :blue, label = "real")
    # plot!(geom.values[nLx,1][ix], imag(P.values)[nLx][ix], color = :red, label  =  "imag")
    # plot!(title = "Straight duct with anechoic termination",
    # xlabel = "x", ylabel = "Pressure")
    # gui()

    true

end
end
using .mstraight_duct_H8_1
mstraight_duct_H8_1.test()

module mmsphere_dipole_1
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsAcoustics
using Test
import LinearAlgebra: norm, lu, cross
function test()

    # println("The interior sphere accelerates in the positive x-direction, generating
    # positive pressure ahead of it, negative pressure behind
    # ")
    rho = 1.21 * phun("kg/m^3")# mass density
    c = 343.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    a_amplitude = 1.0 * phun("mm/s^2")# amplitude of the  acceleration of the sphere
    # omega = 2000*phun("rev/s");      # frequency of the incident wave
    R = 5.0 * phun("mm") # radius of the interior sphere
    Ro = 4 * R # radius of the external sphere
    P_amplitude = R * rho * a_amplitude # pressure amplitude
    nref = 2
    nlayers = 40
    tolerance = Ro / (nlayers) / 100
    frequency = 2000 # Hz
    omega = 2 * pi * frequency
    k = omega / c * [1.0 0.0 0.0]

    # println("""
    #         Spherical domain filled with acoustic medium, no scatterer.
    # Hexahedral mesh.
    #         """)

    t0 = time()

    # Hexahedral mesh
    fens, fes = H8sphere(R, nref)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]],
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)


    # Derive the finite element sets for the boundary
    bfes = meshboundary(fes)
    # Outer spherical boundary
    function dout(xyz)
        return xyz / norm(xyz)
    end

    louter = selectelem(fens, bfes, facing = true, direction = dout)
    outer_fes = subset(bfes, louter)

    fens, fes = H8extrudeQ4(
        fens,
        outer_fes,
        nlayers,
        (xyz, layer) -> (R + layer / nlayers * (Ro - R)) * xyz / norm(xyz),
    )
    connected = findunconnnodes(fens, fes)
    fens, new_numbering = compactnodes(fens, connected)
    fess = renumberconn!(fes, new_numbering)

    geom = NodalField(fens.xyz)
    P = NodalField(fill(zero(FCplxFlt), size(fens.xyz, 1), 1))

    bfes = meshboundary(fes)

    # File  =   "Sphere.vtk"
    # vtkexportmesh(File, bfes.conn, geom.values, FinEtools.MeshExportModule.Q4)
    # @async run(`"paraview.exe" $File`)

    # File  =   "Sphere.vtk"
    # vtkexportmesh (File, fes.conn, fens.xyz, MeshExportModule.H8)
    # @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)

    linner = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0], inflate = tolerance)
    louter = selectelem(fens, bfes, facing = true, direction = dout)

    # println("Pre-processing time elapsed  =  ",time() - t0,"s")

    t1 = time()


    numberdofs!(P)

    material = MatAcoustFluid(bulk, rho)
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material)

    Ma = acousticmass(femm, geom, P)
    Ma = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka = acousticstiffness(femm, geom, P)
    Ka = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]

    abcfemm = FEMMAcoustSurf(IntegDomain(subset(bfes, louter), GaussRule(2, 2)), material)
    D = acousticABC(abcfemm, geom, P)

    # Inner sphere pressure loading
    function dipole(dpdn, xyz, J, label, time = 0.0)
        n = cross(J[:, 1], J[:, 2])
        n = vec(n / norm(n))
        dpdn[1] = -rho * a_amplitude * n[1]
        return dpdn
    end

    fi = ForceIntensity(FCplxFlt, 1, dipole)
    #F  =  F - distribloads(pfemm, nothing, geom, P, fi, 2);
    dipfemm = FEMMAcoustSurf(IntegDomain(subset(bfes, linner), GaussRule(2, 2)), material)
    F = distribloads(dipfemm, geom, P, fi, 2)

    A = (1.0 + 0.0im) * (-omega^2 * Ma + omega * 1.0im * D + Ka)
    K = lu(A) # We fake a complex matrix here
    p = K \ F  #

    scattersysvec!(P, p[:])

    # println(" Minimum/maximum pressure= $(minimum(real(p)))/$(maximum(real(p)))")
    @test abs(-minimum(real(p)) - maximum(real(p))) < 1.0e-9
    @test abs(3.110989923540772e-6 - maximum(real(p))) < 1.0e-9
    # println("Computing time elapsed  =  ",time() - t1,"s")
    # println("Total time elapsed  =  ",time() - t0,"s")

    File = "sphere_dipole_1.vtk"
    vtkexportmesh(
        File,
        connasarray(fes),
        geom.values,
        FinEtools.MeshExportModule.VTK.H8;
        scalars = [("realP", real(P.values))],
    )
    try
        rm(File)
    catch
    end
    # @async run(`"paraview.exe" $File`)

    true
end
end
using .mmsphere_dipole_1
mmsphere_dipole_1.test()

module mstraight_duct_T10_examplem
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsAcoustics
using Test
import Statistics: mean
function test()

    t0 = time()

    rho = 1.21 * phun("kg/m^3")# mass density
    c = 343.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    omega = 54.5901 * phun("rev/s")
    vn0 = -1.0 * phun("m/s")
    Lx = 10.0 * phun("m")# length of the box, millimeters
    Ly = 1.0 * phun("m") # length of the box, millimeters
    n = 20#number of elements along the length

    # println("""
    #
    # Straight duct with anechoic termination.
    # Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 44.
    # Both real and imaginary components of the pressure should have amplitude of
    # rho*c = $(rho*c).
    #
    # Tetrahedral (quadratic) mesh.
    # """)

    fens, fes = T10block(Lx, Ly, Ly, n, 2, 2) # Mesh
    bfes = meshboundary(fes)
    L0 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])
    L10 = selectelem(fens, bfes, facing = true, direction = [+1.0 0.0 0.0])
    nLx = selectnode(fens, box = [0.0 Lx 0.0 0.0 0.0 0.0], inflate = Lx / 1.0e5)

    geom = NodalField(fens.xyz)
    P = NodalField(fill(zero(Complex{Float64}), size(fens.xyz, 1), 1))

    numberdofs!(P)


    material = MatAcoustFluid(bulk, rho)
    femm = FEMMAcoust(IntegDomain(fes, TetRule(4)), material)

    Ma = acousticmass(femm, geom, P)
    Ma = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka = acousticstiffness(femm, geom, P)
    Ka = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]


    E10femm = FEMMAcoustSurf(IntegDomain(subset(bfes, L10), TriRule(3)), material)
    D = acousticABC(E10femm, geom, P)

    E0femm = FEMMBase(IntegDomain(subset(bfes, L0), TriRule(3)))
    fi = ForceIntensity(-1.0im * omega * rho * vn0)
    F = distribloads(E0femm, geom, P, fi, 2)

    p = (-omega^2 * Ma + omega * 1.0im * D + Ka) \ F
    scattersysvec!(P, p[:])

    # println("Pressure amplitude bounds: ")
    # println("  real $(minimum(real(P.values)))/$(maximum(real(P.values)))")
    # println("  imag $(minimum(imag(P.values)))/$(maximum(imag(P.values)))")
    #
    # println("Total time elapsed  =  ",time() - t0,"s")
    @test abs(-414.0986190028914 - minimum(imag(P.values))) < 1.e-6

    File = "straight_duct.vtk"
    scalars = real(P.values)
    vtkexportmesh(
        File,
        connasarray(fes),
        geom.values,
        FinEtools.MeshExportModule.VTK.T10;
        scalars = [("Pressure", scalars)],
    )
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    # plotly()
    # ix = sortperm(geom.values[nLx,1])
    # plot(geom.values[nLx,1][ix], real(P.values)[nLx][ix], color = :blue, label = "real")
    # plot!(geom.values[nLx,1][ix], imag(P.values)[nLx][ix], color = :red, label  =  "imag")
    # plot!(title = "Straight duct with anechoic termination",
    # xlabel = "x", ylabel = "Pressure")
    # gui()

    true
end
end
using .mstraight_duct_T10_examplem
mstraight_duct_T10_examplem.test()

module mmiintegrationmm
using FinEtools
using FinEtoolsAcoustics
using Test
import LinearAlgebra: norm
function test()
    rho = 1000.0 * phun("kg/m^3")# mass density of water
    c = 1.4491e+3 * phun("m/s")# sound speed in water
    bulk = c^2 * rho
    rhos = 2500.0 * phun("kg/m^3")# mass density of the solid sphere
    a_amplitude = 1.0 * phun("mm/s^2")# amplitude of the  acceleration of the sphere
    R = 5.0 * phun("mm") # radius of the interior sphere
    Ro = 3 * R # radius of the external sphere
    # Mesh parameters, geometrical tolerance
    nperR = 6
    nlayers = 10
    tolerance = Ro / (nlayers) / 100
    # Hexahedral mesh of the solid sphere
    fens, fes = H8spheren(R, nperR)
    r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]]
    fens1, fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], renumb = r)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(fens, fes, [0.0, -1.0, 0.0], [0.0, 0.0, 0.0], renumb = r)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(fens, fes, [0.0, 0.0, -1.0], [0.0, 0.0, 0.0], renumb = r)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    # Find the inertial properties of the solid sphere
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    V = integratefunction(femm, geom, (x) -> 1.0)
    # println("V=$(V/phun("mm^3"))")
    Sx = integratefunction(femm, geom, (x) -> x[1])
    # println("Sx=$(Sx/phun("mm^4"))")
    Sy = integratefunction(femm, geom, (x) -> x[2])
    # println("Sy=$(Sy/phun("mm^4"))")
    Sz = integratefunction(femm, geom, (x) -> x[3])
    # println("Sz=$(Sz/phun("mm^4"))")
    CG = vec([Sx Sy Sz] / V)
    # println("CG=$(CG/phun("mm"))")
    I3 = [i == j ? one(FFlt) : zero(FFlt) for i = 1:3, j = 1:3]
    function Iinteg(x)
        r = vec(x) - CG
        (norm(r)^2 * I3 - (r) * (r)')
    end
    Ixx = integratefunction(femm, geom, (x) -> Iinteg(x)[1, 1])
    Ixy = integratefunction(femm, geom, (x) -> Iinteg(x)[1, 2])
    Ixz = integratefunction(femm, geom, (x) -> Iinteg(x)[1, 3])
    Iyy = integratefunction(femm, geom, (x) -> Iinteg(x)[2, 2])
    Iyz = integratefunction(femm, geom, (x) -> Iinteg(x)[2, 3])
    Izz = integratefunction(femm, geom, (x) -> Iinteg(x)[3, 3])
    I = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]
    # println("I=$(I/phun("mm^4"))")
    @test abs(Ixx / phun("mm^4") - 4.9809) < 0.001
    mass = V * rhos
    Inertia = I * rhos
end
end
using .mmiintegrationmm
mmiintegrationmm.test()

module mmtransientsphere
using FinEtools

using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsAcoustics
using Test
import LinearAlgebra: norm, cross
function test()
    # println("The interior sphere accelerates in the alternately in the positive
    # and negative x-direction, generating positive pressure ahead of it, negative
    # pressure behind. Time-dependent simulation.
    # ")
    rho = 1.21 * phun("kg/m^3")# mass density
    c = 343.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    a_amplitude = 1.0 * phun("mm/s^2")# amplitude of the  acceleration of the sphere
    # omega = 2000*phun("rev/s");      # frequency of the incident wave
    R = 50.0 * phun("mm") # radius of the interior sphere
    Ro = 8 * R # radius of the external sphere
    P_amplitude = R * rho * a_amplitude # pressure amplitude
    nref = 2
    nlayers = 40
    tolerance = Ro / (nlayers) / 100
    frequency = 1200.0 # Hz
    omega = 2 * pi * frequency
    dt = 1.0 / frequency / 20
    tfinal = 90 * dt
    nsteps = round(tfinal / dt) + 1

    # Hexahedral mesh
    fens, fes = H8sphere(R, nref)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]],
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    # Derive the finite element sets for the boundary
    bfes = meshboundary(fes)
    # Outer spherical boundary
    function dout(xyz)
        return xyz / norm(xyz)
    end

    louter = selectelem(fens, bfes, facing = true, direction = dout)
    outer_fes = subset(bfes, louter)

    fens, fes = H8extrudeQ4(
        fens,
        outer_fes,
        nlayers,
        (xyz, layer) -> (R + layer / nlayers * (Ro - R)) * xyz / norm(xyz),
    )
    connected = findunconnnodes(fens, fes)
    fens, new_numbering = compactnodes(fens, connected)
    fess = renumberconn!(fes, new_numbering)

    geom = NodalField(fens.xyz)
    P = NodalField(fill(zero(FCplxFlt), size(fens.xyz, 1), 1))

    bfes = meshboundary(fes)

    linner = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0], inflate = tolerance)
    louter = selectelem(fens, bfes, facing = true, direction = dout)

    numberdofs!(P)

    material = MatAcoustFluid(bulk, rho)
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material)

    Ma = acousticmass(femm, geom, P)
    Ma = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka = acousticstiffness(femm, geom, P)
    Ka = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]

    abcfemm = FEMMAcoustSurf(IntegDomain(subset(bfes, louter), GaussRule(2, 2)), material)
    D = acousticABC(abcfemm, geom, P)


    # Inner sphere pressure loading
    function dipole!(dpdn::Vector{FCplxFlt}, xyz, J, feid, qpid, t)
        n = cross(J[:, 1], J[:, 2])
        n = vec(n / norm(n))
        dpdn[1] = -rho * a_amplitude * sin(omega * t) * n[1]
        return dpdn
    end

    dipfemm = FEMMAcoustSurf(IntegDomain(subset(bfes, linner), GaussRule(2, 2)), material)

    # Solve
    P0 = deepcopy(P)
    P0.values[:] .= 0.0 # initially all pressure is zero
    vP0 = gathersysvec(P0)
    vP1 = fill!(similar(vP0), zero(eltype(vP0)))
    vQ0 = fill!(similar(vP0), zero(eltype(vP0)))
    vQ1 = fill!(similar(vP0), zero(eltype(vP0)))

    P1 = deepcopy(P0)

    t = 0.0
    fi = ForceIntensity(
        FCplxFlt,
        1,
        (out, XYZ, tangents, feid, qpid) -> dipole!(out, XYZ, tangents, feid, qpid, t),
    )
    La0 = distribloads(dipfemm, geom, P1, fi, 2)

    A = (2.0 / dt) * Ma + D + (dt / 2.0) * Ka

    step = 0
    while t[] <= tfinal
        step = step + 1
        t = t + dt
        # println("Time $t (  $(step)/$(round(tfinal/dt)+1))")
        fi = ForceIntensity(
            FCplxFlt,
            1,
            (out, XYZ, tangents, feid, qpid) ->
                dipole!(out, XYZ, tangents, feid, qpid, t),
        )
        La1 = distribloads(dipfemm, geom, P1, fi, 2)
        vQ1 =
            A \
            ((2 / dt) * (Ma * vQ0) - D * vQ0 - Ka * (2 * vP0 + (dt / 2) * vQ0) + La0 + La1)
        vP1 = vP0 + (dt / 2) * (vQ0 + vQ1)
        vP0 = vP1
        vQ0 = vQ1
        P1 = scattersysvec!(P1, vec(vP1))
        P0 = P1
        La0 = La1
    end

    @test abs(vP1[1] - (3.0254474104972927e-6 + 0.0im)) < 1.e-10

    # File  =   "sphere_dipole_1.vtk"
    # vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
    #  scalars = [( "realP", real(P1.values))])
    # @async run(`"paraview.exe" $File`)


    true

end
end
using .mmtransientsphere
mmtransientsphere.test()

module mmhhemispheremm
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsAcoustics
using Test
using SparseArrays
import LinearAlgebra: norm, dot, lu, diff, cross
import Statistics: mean
function test()

    # println("Rigid movable hemisphere in  water. Time-dependent simulation.
    # ")
    rho = 1000.0 * phun("kg/m^3")# mass density of water
    c = 1.4491e+3 * phun("m/s")# sound speed in water
    bulk = c^2 * rho
    rhos = 2500.0 * phun("kg/m^3")# mass density of the solid sphere
    a_amplitude = 1.0 * phun("mm/s^2")# amplitude of the  acceleration of the sphere
    R = 5.0 * phun("mm") # radius of the interior sphere
    Ro = 3 * R # radius of the external sphere
    P_amplitude = 1000 * phun("Pa") # pressure amplitude
    frequency = 200.0 # Hz
    omega = 2 * pi * frequency
    wavenumber = omega / c
    angl = 0.0
    front_normal = vec([cos(pi / 180 * angl), sin(pi / 180 * angl), 0])
    wavevector = wavenumber * front_normal
    tshift = 0.0
    uincA = P_amplitude / rho / c / omega # amplitude of the displacement
    # This is the analytical solution of Hickling and Wang
    uamplHW = P_amplitude / (rho * c^2) / (omega / c) * 3 / (2 * rhos / rho + 1)
    # Mesh parameters, geometrical tolerance
    nperR = 8
    nlayers = nperR
    tolerance = Ro / (nlayers) / 100
    # Time step
    dt = 1.0 / frequency / 20
    tfinal = 9 / frequency
    nsteps = Int(round(tfinal / dt) + 1)

    t0 = time()

    # Hexahedral mesh of the solid sphere
    fens, fes = H8spheren(R, nperR)
    r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]]
    fens1, fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], renumb = r)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(fens, fes, [0.0, -1.0, 0.0], [0.0, 0.0, 0.0], renumb = r)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(fens, fes, [0.0, 0.0, -1.0], [0.0, 0.0, 0.0], renumb = r)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    # Retain only half for the hemisphere
    Solidl = selectelem(
        fens,
        fes,
        box = [-Inf, Inf, -Inf, Inf, 0, Inf],
        inflate = (Ro - R) / nlayers / 100,
    )
    # This is the half sphere which is filled with water
    fens1w, fes1w = deepcopy(fens), subset(fes, setdiff(collect(1:count(fes)), Solidl))

    # Find the inertial properties of the solid sphere
    geom = NodalField(fens.xyz)
    femm = FEMMBase(IntegDomain(subset(fes, Solidl), GaussRule(3, 2)))
    V = integratefunction(femm, geom, (x) -> 1.0)
    # println("V=$(V/phun("mm^3"))")
    Sx = integratefunction(femm, geom, (x) -> x[1])
    # println("Sx=$(Sx/phun("mm^4"))")
    Sy = integratefunction(femm, geom, (x) -> x[2])
    # println("Sy=$(Sy/phun("mm^4"))")
    Sz = integratefunction(femm, geom, (x) -> x[3])
    # println("Sz=$(Sz/phun("mm^4"))")
    CG = vec([Sx Sy Sz] / V)
    # println("CG=$(CG/phun("mm"))")
    I3 = [i == j ? one(FFlt) : zero(FFlt) for i = 1:3, j = 1:3]
    function Iinteg(x)
        r = vec(x) - CG
        (norm(r)^2 * I3 - (r) * (r)')
    end
    Ixx = integratefunction(femm, geom, (x) -> Iinteg(x)[1, 1])
    Ixy = integratefunction(femm, geom, (x) -> Iinteg(x)[1, 2])
    Ixz = integratefunction(femm, geom, (x) -> Iinteg(x)[1, 3])
    Iyy = integratefunction(femm, geom, (x) -> Iinteg(x)[2, 2])
    Iyz = integratefunction(femm, geom, (x) -> Iinteg(x)[2, 3])
    Izz = integratefunction(femm, geom, (x) -> Iinteg(x)[3, 3])
    It = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]
    # println("It=$(It/phun("mm^4"))")
    mass = V * rhos
    Inertia = It * rhos

    # The boundary surface of the solid sphere will now be extruded to form  the
    # layers of fluid around it
    # Outer spherical boundary
    function dout(xyz)
        return xyz / norm(xyz)
    end
    bfes = meshboundary(fes)
    louter = selectelem(fens, bfes, facing = true, direction = dout)
    outer_fes = subset(bfes, louter)
    fens, fes = H8extrudeQ4(
        fens,
        outer_fes,
        nlayers,
        (xyz, layer) -> (R + layer / nlayers * (Ro - R)) * xyz / norm(xyz),
    )
    fens, newfes1, fes2 = mergemeshes(fens, fes, fens1w, fes1w, tolerance)
    fes = cat(newfes1, fes2)
    connected = findunconnnodes(fens, fes)
    fens, new_numbering = compactnodes(fens, connected)
    fes = renumberconn!(fes, new_numbering)


    # File  =   "Sphere1.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)

    geom = NodalField(fens.xyz)
    P = NodalField(fill(zero(FFlt), size(fens.xyz, 1), 1))

    bfes = meshboundary(fes)

    linner = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0], inflate = tolerance)
    louter = selectelem(fens, bfes, facing = true, direction = dout)

    # println("Pre-processing time elapsed  =  ",time() - t0,"s")

    t1 = time()

    numberdofs!(P)

    material = MatAcoustFluid(bulk, rho)
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material)

    Ma = acousticmass(femm, geom, P)
    Ma = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka = acousticstiffness(femm, geom, P)
    Ka = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]

    abcfemm = FEMMAcoustSurf(IntegDomain(subset(bfes, louter), GaussRule(2, 2)), material)
    D = acousticABC(abcfemm, geom, P)

    # ABC surface pressure loading
    function abcp(dpdn, xyz, J, t)
        n = cross(J[:, 1], J[:, 2])
        n = vec(n / norm(n))
        arg = (-dot(vec(xyz), vec(wavevector)) + omega * (t + tshift))
        dpdn[1] = P_amplitude * cos(arg) * (-dot(vec(n), vec(wavevector)))
        return dpdn
    end


    targetfemm =
        FEMMAcoustSurf(IntegDomain(subset(bfes, linner), GaussRule(2, 2)), material)

    ForceF = GeneralField([zero(FFlt) for i = 1:3, j = 1:1])
    numberdofs!(ForceF)
    TorqueF = GeneralField([zero(FFlt) for i = 1:3, j = 1:1])
    numberdofs!(TorqueF)

    GF = pressure2resultantforce(targetfemm, geom, P, ForceF)
    GT = pressure2resultanttorque(targetfemm, geom, P, TorqueF, CG)

    H =
        transpose(GF) * (rho / mass) * GF +
        transpose(GT) * (sparse(rho * inv(Inertia))) * GT

    Ctild = Ka + H

    # Solve
    P0 = deepcopy(P)
    P0.values[:] .= 0.0 # initially all pressure is zero
    vP0 = gathersysvec(P0)
    vP1 = fill!(similar(vP0), zero(eltype(vP0)))
    vQ0 = fill!(similar(vP0), zero(eltype(vP0)))
    vQ1 = fill!(similar(vP0), zero(eltype(vP0)))
    t = 0.0
    P1 = deepcopy(P0)
    tTa1 = [zero(FFlt) for i = 1:3]
    tAa1 = [zero(FFlt) for i = 1:3]
    tTa0 = [zero(FFlt) for i = 1:3]
    tAa0 = [zero(FFlt) for i = 1:3]
    tTv1 = [zero(FFlt) for i = 1:3]
    tAv1 = [zero(FFlt) for i = 1:3]
    tTv0 = [zero(FFlt) for i = 1:3]
    tAv0 = [zero(FFlt) for i = 1:3]
    tTu1 = [zero(FFlt) for i = 1:3]
    tAu1 = [zero(FFlt) for i = 1:3]
    tTu0 = [zero(FFlt) for i = 1:3]
    tAu0 = [zero(FFlt) for i = 1:3]
    tTa_store = fill(zero(FFlt), 3, nsteps + 1)
    tAa_store = fill!(similar(tTa_store), zero(eltype(tTa_store)))
    tTv_store = fill!(similar(tTa_store), zero(eltype(tTa_store)))
    tAv_store = fill!(similar(tTa_store), zero(eltype(tTa_store)))
    tTu_store = fill!(similar(tTa_store), zero(eltype(tTa_store)))
    tAu_store = fill!(similar(tTa_store), zero(eltype(tTa_store)))
    t_store = fill(zero(FFlt), 1, nsteps + 1)

    pinc = deepcopy(P)
    pincdd = deepcopy(P)
    pincv = fill(zero(FFlt), nalldofs(P))
    pincddv = deepcopy(pincv)
    p1v = deepcopy(pincv)
    p1 = deepcopy(P)

    function recalculate_incident!(t, pinc, pincdd)
        for j = 1:count(fens)
            arg = (-dot(vec(fens.xyz[j, :]), vec(wavevector)) + omega * (t + tshift))# wave along the wavevector
            pinc.values[j] = P_amplitude * sin(arg)
            pincdd.values[j] = -omega^2 * P_amplitude * sin(arg)
        end
        return pinc, pincdd
    end

    pinc, pincdd = recalculate_incident!(t, pinc, pincdd)
    gathersysvec!(pincdd, pincddv)
    gathersysvec!(pinc, pincv)
    L0 = -Ma * pincddv[freedofs(pincdd)] - Ctild * pincv[freedofs(pincdd)]
    fi = ForceIntensity(FFlt, 1, (dpdn, xyz, J, feid, qpid) -> abcp(dpdn, xyz, J, t))
    La0 = distribloads(abcfemm, geom, P1, fi, 2)

    A = lu((2.0 / dt) * Ma + D + (dt / 2.0) * Ctild)

    step = 1
    while step <= nsteps
        # println("Time $t ($(step)/$(Int(round(tfinal/dt))+1))")
        # Store for output
        tTa_store[:, step] = tTa1
        tAa_store[:, step] = tAa1
        tTv_store[:, step] = tTv1
        tAv_store[:, step] = tAv1
        tTu_store[:, step] = tTu1
        tAu_store[:, step] = tAu1
        t_store[1, step] = t
        # New loads and update
        step = step + 1
        t = t + dt
        pinc, pincdd = recalculate_incident!(t, pinc, pincdd)
        gathersysvec!(pincdd, pincddv)
        gathersysvec!(pinc, pincv)
        L1 = -Ma * pincddv[freedofs(pincdd)] - Ctild * pincv[freedofs(pincdd)]
        fi = ForceIntensity(FFlt, 1, (dpdn, xyz, J, feid, qpid) -> abcp(dpdn, xyz, J, t))
        La1 = distribloads(abcfemm, geom, P1, fi, 2)
        vQ1 =
            A \ (
                (2 / dt) * (Ma * vQ0) - D * vQ0 - Ctild * (2 * vP0 + (dt / 2) * vQ0) +
                L0 +
                L1 +
                La0 +
                La1
            )
        vP1 = vP0 + (dt / 2) * (vQ0 + vQ1)
        P1 = scattersysvec!(P1, vec(vP1))
        p1.values[:] = P1.values[:] + pinc.values[:]
        p1v = gathersysvec!(p1, p1v)
        Fresultant = GF * p1v
        Tresultant = GT * p1v
        # println("Fresultant = $Fresultant")
        tTa1 = Fresultant / mass
        tAa1 = Inertia \ Tresultant
        # Update  the velocity and position of the rigid body
        tTv1 = tTv0 + dt / 2 * (tTa0 + tTa1)
        tTu1 = tTu0 + dt / 2 * (tTv0 + tTv1)
        tAv1 = tAv1 + dt / 2 * (tAa0 + tAa1)
        tAu1 = tAu0 + dt / 2 * (tAv0 + tAv1)
        # show(tTa1)
        # Reset for the next time step
        copyto!(vP0, vP1)
        copyto!(vQ0, vQ1)
        P0 = deepcopy(P1)
        copyto!(L0, L1)
        La0, La1 = La1, La0
        copyto!(tTa0, tTa1)
        copyto!(tAa0, tAa1)
        copyto!(tTv0, tTv1)
        copyto!(tAv0, tAv1)
        copyto!(tTu0, tTu1)
        copyto!(tAu0, tAu1)
    end
    tTa_store[:, step] = tTa1
    tAa_store[:, step] = tAa1
    tTv_store[:, step] = tTv1
    tAv_store[:, step] = tAv1
    tTu_store[:, step] = tTu1
    tAu_store[:, step] = tAu1
    t_store[1, step] = t

    function remove_bias(it, iy)
        local t = deepcopy(reshape(it, length(it), 1))
        local y = deepcopy(reshape(iy, length(iy), 1))
        t = t / mean(t)
        t = reshape(t, length(t), 1)
        y = reshape(y, length(y), 1)
        local B = [t [1.0 for m in t]]
        local p = (B' * B) \ (B' * y)
        return y .- p[1] * t .- p[2]
    end

    ux = remove_bias(t_store, tTu_store[1, :])

    function find_peaks(t, y)
        p = Float64[]
        dy = diff(y, dims = 1)
        for j3 = 1:length(dy)-1
            if (dy[j3] * dy[j3+1] < 0)
                push!(p, y[j3+1])
            end
        end
        return p
    end
    function amplitude(t, y)
        p = find_peaks(t, y)
        A = mean(abs.(p))
        return A
    end
    # println("Normalized amplitude  = $(amplitude(t_store,ux)/uamplHW)")
    @test abs(amplitude(t_store, ux) / uamplHW - 0.9361240886218495) < 1.0e-5
    # using Plots
    # plotly()
    # plot(transpose(t_store), transpose(tTu_store)/phun("mm"),
    #  xlabel = "t", ylabel = "Displacement", label = "Displacements")
    # plot(transpose(t_store), transpose(tAu_store),
    #  xlabel = "t", ylabel = "Displacement", label = "Rotations")
    # # plot(transpose(t_store), (ux)/uamplHW,
    # #  xlabel = "t", ylabel = "Normalized displacement", label = "Rectified")
    # gui()
    #
    # #
    # File  =   "hemisphere_final_td1.vtk"
    # vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
    #  scalars = [( "realP", real(P1.values))])
    # @async run(`"paraview.exe" $File`)

end
end
using .mmhhemispheremm
mmhhemispheremm.test()

module mmbbaffledmm
using FinEtools
using FinEtoolsAcoustics
using FinEtoolsAcoustics.AlgoAcoustModule
using Test
import LinearAlgebra: norm
function test()
    rho = 1.21 * phun("kg/m^3")# mass density
    c = 343.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    omega = 7500 * phun("rev/s")      # frequency of the piston
    a_piston = -1.0 * phun("mm/s")     # amplitude of the piston acceleration
    R = 50.0 * phun("mm")# radius of the piston
    Ro = 150.0 * phun("mm") # radius of the enclosure
    nref = 3#number of refinements of the sphere around the piston
    nlayers = 25                     # number of layers of elements surrounding the piston
    tolerance = R / (2^nref) / 100

    # println("""
    #
    # Baffled piston in a half-sphere domain with ABC.
    #
    # Hexahedral mesh. Algorithm version.
    # """)

    t0 = time()

    # Hexahedral mesh
    fens, fes = H8sphere(R, nref)
    bfes = meshboundary(fes)
    # File  =   "baffledabc_boundary.vtk"
    # vtkexportmesh(File, bfes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
    #  @async run(`"paraview.exe" $File`)

    l = selectelem(fens, bfes, facing = true, direction = [1.0 1.0 1.0], dotmin = 0.001)
    ex(xyz, layer) = (R + layer / nlayers * (Ro - R)) * xyz / norm(xyz)
    fens1, fes1 = H8extrudeQ4(fens, subset(bfes, l), nlayers, ex)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    # Piston surface mesh
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])
    l2 = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0], inflate = tolerance)
    piston_fes = subset(bfes, intersect(l1, l2))

    # Outer spherical boundary
    louter =
        selectelem(fens, bfes, facing = true, direction = [1.0 1.0 1.0], dotmin = 0.001)
    outer_fes = subset(bfes, louter)

    # println("Pre-processing time elapsed  =  ",time() - t0,"s")

    t1 = time()

    material = MatAcoustFluid(bulk, rho)
    # Region of the fluid
    region1 = FDataDict("femm" => FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material))

    # Surface for the ABC
    abc1 = FDataDict(
        "femm" => FEMMAcoustSurf(IntegDomain(outer_fes, GaussRule(2, 2)), material),
    )

    # Surface of the piston
    flux1 = FDataDict(
        "femm" => FEMMAcoustSurf(IntegDomain(piston_fes, GaussRule(2, 2)), material),
        "normal_flux" => -1.0im * rho * a_piston,
    )

    # Make model data
    modeldata = FDataDict(
        "fens" => fens,
        "omega" => omega,
        "regions" => [region1],
        "flux_bcs" => [flux1],
        "ABCs" => [abc1],
    )

    # Call the solver
    modeldata = AlgoAcoustModule.steadystate(modeldata)

    # println("Computing time elapsed  =  ",time() - t1,"s")
    # println("Total time elapsed  =  ",time() - t0,"s")

    geom = modeldata["geom"]
    P = modeldata["P"]

    @test abs(P.values[1] - (1.2061599990626182e-6 + 4.355914465356351e-6im)) < 1.e-10
    # File  =   "baffledabc.vtk"
    # vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
    # scalars = [("absP", abs.(P.values))])
    # @async run(`"paraview.exe" $File`)

    # using Winston
    # pl  =  FramedPlot(title = "Matrix",xlabel = "x",ylabel = "Re P, Im P")
    # setattr(plrame, draw_grid = true)
    # add(pl, Curve([1:length(C[:])],vec(C[:]), color = "blue"))

    # # pl = plot(geom.values[nLx,1][ix],scalars[nLx][ix])
    # # xlabel("x")
    # # ylabel("Pressure")
    # display(pl)

    true
end
end
using .mmbbaffledmm
mmbbaffledmm.test()

module mmtransientmm1mm
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsAcoustics
using Test
import LinearAlgebra: norm
function test()
    rho = 1.21 * phun("kg/m^3")# mass density
    c = 343.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    frequency = 2000 # Hz
    omega = 2 * pi * frequency
    P_piston = 1.0e-3 * phun("MPa") # amplitude of the piston pressure
    R = 50.0 * phun("mm")# radius of the piston
    Ro = 400.0 * phun("mm") # radius of the enclosure
    nref = 2#number of refinements of the sphere around the piston
    nlayers = 48                     # number of layers of elements surrounding the piston
    tolerance = R / (2^nref) / 10
    dt = 1.0 / frequency / 20
    tfinal = 50 * dt
    nsteps = round(tfinal / dt) + 1

    # println("""
    #
    # Baffled piston in a half-sphere domain with hard boundary condition.
    #
    # Hexahedral mesh. Algorithm version.
    # """)

    t0 = time()

    # Hexahedral mesh
    fens, fes = H8sphere(R, nref)
    bfes = meshboundary(fes)
    # File  =   "baffledabc_boundary.vtk"
    # vtkexportmesh(File, bfes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
    #  @async run(`"paraview.exe" $File`)

    l = selectelem(fens, bfes, facing = true, direction = [1.0 1.0 1.0], dotmin = 0.001)
    ex(xyz, layer) = (R + layer / nlayers * (Ro - R)) * xyz / norm(xyz)
    fens1, fes1 = H8extrudeQ4(fens, subset(bfes, l), nlayers, ex)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    # Piston surface mesh
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])
    l2 = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0], inflate = tolerance)
    piston_fes = subset(bfes, intersect(l1, l2))
    # File  =   "piston_fes.vtk"
    # vtkexportmesh(File, piston_fes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
    #  @async run(`"paraview.exe" $File`)

    # Outer spherical boundary
    louter =
        selectelem(fens, bfes, facing = true, direction = [1.0 1.0 1.0], dotmin = 0.001)
    outer_fes = subset(bfes, louter)
    # File  =   "outer_fes.vtk"
    # vtkexportmesh(File, outer_fes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
    #  @async run(`"paraview.exe" $File`)

    # println("Pre-processing time elapsed  =  ",time() - t0,"s")

    t1 = time()

    material = MatAcoustFluid(bulk, rho)
    # Region of the fluid
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material)

    geom = NodalField(fens.xyz)
    P = NodalField(fill(zero(FFlt), nnodes(geom), 1))
    piston_fenids = connectednodes(piston_fes)
    setebc!(P, piston_fenids, true, 1, 0.0)
    applyebc!(P)
    numberdofs!(P)

    Ma = acousticmass(femm, geom, P)
    Ka = acousticstiffness(femm, geom, P)

    Ma_ff, Ma_fd = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[(:ff, :fd)]
    Ka_ff, Ka_fd = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[(:ff, :fd)]

    F = zeros(ComplexF64, nalldofs(P))
    F_f, F_d = vector_blocked(F, nfreedofs(P))[(:f, :d)]

    A = (2.0 / dt) * Ma_ff + (dt / 2.0) * Ka_ff

    P0 = deepcopy(P)
    Pdd0 = deepcopy(P)
    P1 = deepcopy(P)
    Pdd1 = deepcopy(P)
    TMPF = deepcopy(P)
    P = nothing # we don't need this field anymore
    vP0 = fill(zero(FFlt), nfreedofs(P0))
    vP1 = deepcopy(vP0)
    vPd0 = deepcopy(vP0)
    vPd1 = deepcopy(vP0)

    t = 0.0
    P0.values[piston_fenids, 1] .= P_piston * sin(omega * t)
    Pdd0.values[piston_fenids, 1] .= P_piston * (-omega^2) * sin(omega * t)
    vP0 = gathersysvec!(P0, vP0)[freedofs(P0)]

    nh = selectnode(fens, nearestto = [R + Ro / 2, 0.0, 0.0])
    Pnh = [P1.values[nh, 1][1]]
    ts = [t]
    step = 0
    while t <= tfinal
        step = step + 1
        t = t + dt
        P1.values[piston_fenids, 1] .= P_piston * sin(omega * t)
        Pdd1.values[piston_fenids, 1] .= P_piston * (-omega^2) * sin(omega * t)
        TMPF.values = P0.values + P1.values
        gathersysvec!(TMPF, F_d, :d)
        F_f .= -Ka_fd * F_d
        TMPF.values = Pdd0.values + Pdd1.values
        gathersysvec!(TMPF, F_d, :d)
        F_f = F_f - Ma_fd * F_d
        # println("$(norm(F))")
        vPd1 = A \ ((2 / dt) * (Ma_ff * vPd0) - Ka_ff * (2 * vP0 + (dt / 2) * vPd0) + F_f)
        vP1 = vP0 + (dt / 2) * (vPd0 + vPd1)
        scattersysvec!(P1, vP1) # store current pressure
        # Swap variables for the next step
        copyto!(vP0, vP1)
        copyto!(vPd0, vPd1)
        copyto!(P0, P1)
        copyto!(Pdd0, Pdd1)
        # Graphics output
        # File  =   "baffled_piston-$(step).vtk"
        # vtkexportmesh(File, fes.conn, fens.xyz, FinEtools.MeshExportModule.H8;
        #     scalars = [(  "P", P1.values)])
        #  @async run(`"paraview.exe" $File`)
        push!(Pnh, P1.values[nh, 1][1])
        push!(ts, t)
        # println("step = $( step )")
    end

    # using Plots
    # plotly()
    # plot(vec(ts), vec(Pnh))
    # gui()
    # println("$(Pnh[end])")
    # println("$(length(Pnh))")
    # println("$(Pnh[end])")
    @test abs(Pnh[end] - 221.11319820621947) < 1.0e-3
end
end
using .mmtransientmm1mm
mmtransientmm1mm.test()

module mmAnnularmm
using FinEtools
using FinEtoolsAcoustics
using FinEtoolsAcoustics.AlgoAcoustModule
using Test
function test()


    # println("""
    # Annular region, pressure BC + rigid wall.
    # This version uses the FinEtools algorithm module.
    # Version: 08/21/2017
    # """)

    t0 = time()

    rho = 1001 * phun("kg/m^3")# mass density
    c = 1500.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    omega = 2000 * phun("rev/s")      # frequency of the piston
    rin = 1.0 * phun("m")#internal radius

    rex = 2.0 * phun("m") #external radius
    nr = 20
    nc = 120
    Angle = 2 * pi
    thickness = 1.0 * phun("m/s")
    tolerance = min(rin / nr, rin / nc / 2 / pi) / 10000

    fens, fes = Q4annulus(rin, rex, nr, nc, Angle)
    fens, fes = mergenodes(fens, fes, tolerance)
    edge_fes = meshboundary(fes)

    # The pressure boundary condition
    l1 = selectelem(fens, edge_fes, box = [-1.1 * rex -0.9 * rex -0.5 * rex 0.5 * rex])
    ebc1 = FDataDict(
        "node_list" => connectednodes(subset(edge_fes, l1)),
        "pressure" => x -> cos(2 * pi * x[2] / rin) + 1im * sin(2 * pi * x[2] / rin),
    ) # entering the domain

    material = MatAcoustFluid(bulk, rho)
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(2, 2)), material)
    region1 = FDataDict("femm" => femm)

    # Make model data
    modeldata = FDataDict(
        "fens" => fens,
        "omega" => omega,
        "regions" => [region1],
        "essential_bcs" => [ebc1],
    )

    # Call the solver
    modeldata = AlgoAcoustModule.steadystate(modeldata)
    geom = modeldata["geom"]
    P = modeldata["P"]
    # println("Minimum/maximum pressure, real= $(minimum(real(P.values)))/$(maximum(real(P.values))))")
    # println("Minimum/maximum pressure, imag= $(minimum(imag(P.values)))/$(maximum(imag(P.values))))")
    @test abs(minimum(real(P.values)) - -26.02026924534437) < 1.0e-5
    @test abs(maximum(real(P.values)) - 32.15639028990587) < 1.0e-5
    @test abs(minimum(imag(P.values)) - -3.083023162317023) < 1.0e-5
    @test abs(maximum(imag(P.values)) - 3.083023162317045) < 1.0e-5
    # println("Total time elapsed = ",time() - t0,"s")

    # # Postprocessing
    # File = "acou_annulusmod.vtk"
    # vtkexportmesh(File, fes.conn, geom.values,
    # FinEtools.MeshExportModule.Q4; scalars=[("Pre", real(P.values)), ("Pim", imag(P.values))])
    # @async run(`"paraview.exe" $File`)

end
end
using .mmAnnularmm
mmAnnularmm.test()

module mmbbaffledmAlgo
using FinEtools
using FinEtoolsAcoustics
using Test
import LinearAlgebra: norm
function test()
    rho = 1.21 * phun("kg/m^3")# mass density
    c = 343.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    omega = 7500 * phun("rev/s")      # frequency of the piston
    a_piston = -1.0 * phun("mm/s")     # amplitude of the piston acceleration
    R = 50.0 * phun("mm")# radius of the piston
    Ro = 150.0 * phun("mm") # radius of the enclosure
    nref = 3#number of refinements of the sphere around the piston
    nlayers = 25                     # number of layers of elements surrounding the piston
    tolerance = R / (2^nref) / 100

    # println("""
    #
    # Baffled piston in a half-sphere domain with ABC.
    #
    # Hexahedral mesh. Algorithm version.
    # """)

    t0 = time()

    # Hexahedral mesh
    fens, fes = H8sphere(R, nref)
    bfes = meshboundary(fes)
    # File  =   "baffledabc_boundary.vtk"
    # vtkexportmesh(File, bfes.conn, fens.xyz, FinEtools.MeshExportModule.Q4)
    #  @async run(`"paraview.exe" $File`)

    l = selectelem(fens, bfes, facing = true, direction = [1.0 1.0 1.0], dotmin = 0.001)
    ex(xyz, layer) = (R + layer / nlayers * (Ro - R)) * xyz / norm(xyz)
    fens1, fes1 = H8extrudeQ4(fens, subset(bfes, l), nlayers, ex)
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    # Piston surface mesh
    bfes = meshboundary(fes)
    l1 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])
    l2 = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0], inflate = tolerance)
    piston_fes = subset(bfes, intersect(l1, l2))

    # Outer spherical boundary
    louter =
        selectelem(fens, bfes, facing = true, direction = [1.0 1.0 1.0], dotmin = 0.001)
    outer_fes = subset(bfes, louter)

    # println("Pre-processing time elapsed  =  ",time() - t0,"s")

    t1 = time()

    material = MatAcoustFluid(bulk, rho)
    # Region of the fluid
    region1 = FDataDict("femm" => FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material))

    # Surface for the ABC
    abc1 = FDataDict(
        "femm" => FEMMAcoustSurf(IntegDomain(outer_fes, GaussRule(2, 2)), material),
    )

    # Surface of the piston
    flux1 = FDataDict(
        "femm" => FEMMAcoustSurf(IntegDomain(piston_fes, GaussRule(2, 2)), material),
        "normal_flux" =>
            (forceout, XYZ, tangents, feid, qpid) ->
                (forceout[1] = -1.0im * rho * a_piston; forceout),
    )

    # Make model data
    modeldata = FDataDict(
        "fens" => fens,
        "omega" => omega,
        "regions" => [region1],
        "flux_bcs" => [flux1],
        "ABCs" => [abc1],
    )

    # Call the solver
    modeldata = FinEtoolsAcoustics.AlgoAcoustModule.steadystate(modeldata)

    # println("Computing time elapsed  =  ",time() - t1,"s")
    # println("Total time elapsed  =  ",time() - t0,"s")

    geom = modeldata["geom"]
    P = modeldata["P"]

    @test abs(P.values[1] - (4.355914465396856e-6im + 1.2061599990585114e-6)) < 1.e-10
    # File  =   "baffledabc.vtk"
    # vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
    # scalars = [("absP", abs.(P.values))])
    # @async run(`"paraview.exe" $File`)

    # using Winston
    # pl  =  FramedPlot(title = "Matrix",xlabel = "x",ylabel = "Re P, Im P")
    # setattr(pl.frame, draw_grid = true)
    # add(pl, Curve([1:length(C[:])],vec(C[:]), color = "blue"))

    # # pl = plot(geom.values[nLx,1][ix],scalars[nLx][ix])
    # # xlabel("x")
    # # ylabel("Pressure")
    # display(pl)

    true
end
end
using .mmbbaffledmAlgo
mmbbaffledmAlgo.test()


module mmAnnularmAlgo
using FinEtools
using FinEtoolsAcoustics
using FinEtoolsAcoustics.AlgoAcoustModule
using Test
function test()


    # println("""
    # Annular region, pressure BC + rigid wall.
    # This version uses the FinEtools algorithm module.
    # Version: 08/21/2017
    # """)

    t0 = time()

    rho = 1001 * phun("kg/m^3")# mass density
    c = 1500.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    omega = 2000 * phun("rev/s")      # frequency of the piston
    rin = 1.0 * phun("m")#internal radius

    rex = 2.0 * phun("m") #external radius
    nr = 20
    nc = 120
    Angle = 2 * pi
    thickness = 1.0 * phun("m/s")
    tolerance = min(rin / nr, rin / nc / 2 / pi) / 10000

    fens, fes = Q4annulus(rin, rex, nr, nc, Angle)
    fens, fes = mergenodes(fens, fes, tolerance)
    edge_fes = meshboundary(fes)

    # The pressure boundary condition
    l1 = selectelem(fens, edge_fes, box = [-1.1 * rex -0.9 * rex -0.5 * rex 0.5 * rex])
    ebc1 = FDataDict(
        "node_list" => connectednodes(subset(edge_fes, l1)),
        "pressure" => 1.0 + 0.0im,
    ) # entering the domain

    material = MatAcoustFluid(bulk, rho)
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(2, 2)), material)
    region1 = FDataDict("femm" => femm)

    # Make model data
    modeldata = FDataDict(
        "fens" => fens,
        "omega" => omega,
        "regions" => [region1],
        "essential_bcs" => [ebc1],
    )

    # Call the solver
    modeldata = AlgoAcoustModule.steadystate(modeldata)
    geom = modeldata["geom"]
    P = modeldata["P"]
    # println("Minimum/maximum pressure, real= $(minimum(real(P.values)))/$(maximum(real(P.values))))")
    # println("Minimum/maximum pressure, imag= $(minimum(imag(P.values)))/$(maximum(imag(P.values))))")
    @test abs(minimum(real(P.values)) - -97.34672558165316) < 1.0e-5
    @test abs(maximum(real(P.values)) - 75.50540650112683) < 1.0e-5
    @test abs(minimum(imag(P.values)) - 0.0) < 1.0e-5
    @test abs(maximum(imag(P.values)) - 0.0) < 1.0e-5
    # println("Total time elapsed = ",time() - t0,"s")

    # # Postprocessing
    # File = "acou_annulusmod.vtk"
    # vtkexportmesh(File, fes.conn, geom.values,
    # FinEtools.MeshExportModule.Q4; scalars=[("Pre", real(P.values)), ("Pim", imag(P.values))])
    # @async run(`"paraview.exe" $File`)

end
end
using .mmAnnularmAlgo
mmAnnularmAlgo.test()

module mmgradient_1
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsAcoustics
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using Test
import LinearAlgebra: cholesky

function h20_mesh()
    Lx = 0.414 * phun("m") # dimension
    Ly = 0.314 * phun("m") # dimension
    Lz = 0.360 * phun("m") # dimension
    n = 3
    fens, fes = H8block(Lx, Ly, Lz, n, n, n) # Mesh
    fens, fes = H8toH20(fens, fes)
    fens, fes
end

function h8_mesh()
    Lx = 0.414 * phun("m") # dimension
    Ly = 0.314 * phun("m") # dimension
    Lz = 0.360 * phun("m") # dimension
    n = 10
    fens, fes = H8block(Lx, Ly, Lz, n, n, n) # Mesh
    fens, fes
end

const meshfun = Dict("h20" => h20_mesh, "h8" => h8_mesh)

function test(name)
    rho = 1.21 * phun("kg/m^3") # mass density
    c = 343.0 * phun("m/s")  # sound speed
    bulk = c^2 * rho
    alpha, beta, gamma = 0.13, -0.3, 0.1333
    pfun(x) = alpha * x[1] + beta * x[2] + gamma * x[3]

    fens, fes = meshfun[name]()

    geom = NodalField(fens.xyz)
    P = NodalField(reshape([pfun(fens.xyz[i, :]) for i in eachindex(fens)], count(fens), 1))
    numberdofs!(P)

    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 3)), MatAcoustFluid(bulk, rho))

    gradP = reshape([alpha, beta, gamma], 1, 3)
    function inspector(idat, i, conn, xe, out, loc)
        erro = maximum(abs.(gradP .- out))
        return erro
    end

    idat = inspectintegpoints(femm, geom, NodalField([1.0]), P, collect(1:count(fes)), inspector, 0.0, :gradient)

    @test idat  < 1.0e-9 / count(fes)
    true
end #

test("h8")
test("h20")

end # module
nothing

module mmgradient_2
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsAcoustics
using FinEtools.MeshExportModule.VTK: vtkexportmesh, T3, vtkexportvectors
using Test
import LinearAlgebra: cholesky

function h20_mesh()
    Lx = 0.414 * phun("m") # dimension
    Ly = 0.314 * phun("m") # dimension
    Lz = 0.360 * phun("m") # dimension
    n = 3
    fens, fes = H8block(Lx, Ly, Lz, n, n, n) # Mesh
    fens, fes = H8toH20(fens, fes)
    fens, fes
end

function h8_mesh()
    Lx = 0.414 * phun("m") # dimension
    Ly = 0.314 * phun("m") # dimension
    Lz = 0.360 * phun("m") # dimension
    n = 10
    fens, fes = H8block(Lx, Ly, Lz, n, n, n) # Mesh
    fens, fes
end

const meshfun = Dict("h20" => h20_mesh, "h8" => h8_mesh)

function test(name)
    rho = 1.21 * phun("kg/m^3") # mass density
    c = 343.0 * phun("m/s")  # sound speed
    bulk = c^2 * rho
    alpha, beta, gamma = 0.13, -0.3, 0.1333
    pfun(x) = alpha * x[1] + beta * x[2] + gamma * x[3]
    gradP = reshape([alpha, beta, gamma], 1, 3)

    fens, fes = meshfun[name]()

    geom = NodalField(fens.xyz)
    P = NodalField(reshape([pfun(fens.xyz[i, :]) for i in eachindex(fens)], count(fens), 1))
    numberdofs!(P)

    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 3)), MatAcoustFluid(bulk, rho))

    qplocs = []
    qpgradients = []
    function inspector(idat, i, conn, xe, out, loc)
      qplocs, qpgradients = idat
      push!(qplocs, copy(loc))
      push!(qpgradients, copy(out))
      return (qplocs, qpgradients)
    end

    idat = inspectintegpoints(femm, geom, NodalField([1.0]), P, collect(1:count(fes)), inspector, (qplocs, qpgradients), :gradient)
    qplocs, qpgradients = idat

    # Test export of vectors of gradient
    File =  "mmgradient_2-$(name)-vectors.vtk"
    vtkexportvectors(File, qplocs, [("grad", qpgradients)])
    try rm(File); catch end
    File =  "mmgradient_2-$(name).vtk"
    vtkexportmesh(File, fens, fes)
    try rm(File); catch end

    true
end #

test("h8")
test("h20")

end # module
nothing

