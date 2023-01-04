module Fahy_examples
using FinEtools
using FinEtoolsAcoustics
using PlotlyLight
import Arpack: eigs

function fahy_H20_example()
    println("""
    Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
    Frank J. Fahy, Paolo Gardonio, page 483.
    
    Hexahedral mesh.
    """)
    
    t0 = time()
    
    rho=1.21*1e-9;# mass density
    c =343.0*1000;# millimeters per second
    bulk= c^2*rho;
    L=500.0;# length of the box, millimeters
    A=200.0; # cross-sectional area of the box
    n=40;#
    neigvs=8;
    OmegaShift=10.0;
    
    fens,fes = H8block(L,sqrt(A),sqrt(A),n,1,1); # Mesh
    fens,fes = H8toH20(fens,fes)
    
    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz,1),1))
    
    numberdofs!(P)
    
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 3)), MatAcoustFluid(bulk, rho))
    
    
    S = acousticstiffness(femm, geom, P);
    C = acousticmass(femm, geom, P);
    
    d,v,nev,nconv = eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    v = real.(v)
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    
    println("Total time elapsed = ",time() - t0,"s")
    
    File =  "fahy_H20.vtk"
    en = 5;
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.VTK.H20;
    scalars=[("Pressure_mode_$en", v[:,en])])
    # @async run(`"paraview.exe" $File`)
    println("Done")
    true
    
end # fahy_H20_example


function fahy_H27_example()
    println("""
    Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
    Frank J. Fahy, Paolo Gardonio, page 483.
    
    Hexahedral mesh. 27-node elements.
    """)
    
    t0 = time()
    
    rho=1.21*1e-9;# mass density
    c =343.0*1000;# millimeters per second
    bulk= c^2*rho;
    L=500.0;# length of the box, millimeters
    A=200.0; # cross-sectional area of the box
    n=40;#
    neigvs=8;
    OmegaShift=10.0;
    
    fens,fes = H27block(L,sqrt(A),sqrt(A),n,1,1); # Mesh
    
    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz,1),1))
    
    numberdofs!(P)
    
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 3)), MatAcoustFluid(bulk, rho))
    
    
    S = acousticstiffness(femm, geom, P);
    C = acousticmass(femm, geom, P);
    
    d,v,nev,nconv = eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    v = real.(v)
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    
    println("Total time elapsed = ",time() - t0,"s")
    
    
    println("Done")
    true
    
end # fahy_H27_example


function fahy_H8_example()
    println("""
    Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
    Frank J. Fahy, Paolo Gardonio, page 483.
    
    Hexahedral mesh.
    """)
    
    t0 = time()
    
    rho=1.21*1e-9;# mass density
    c =343.0*1000;# millimeters per second
    bulk= c^2*rho;
    L=500.0;# length of the box, millimeters
    A=200.0; # cross-sectional area of the box
    n=40;#
    neigvs=8;
    OmegaShift=10.0;
    
    fens,fes = H8block(L,sqrt(A),sqrt(A),n,1,1); # Mesh
    
    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz,1),1))
    
    numberdofs!(P)
    
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), MatAcoustFluid(bulk, rho))
    
    
    S = acousticstiffness(femm, geom, P);
    C = acousticmass(femm, geom, P);
    
    d,v,nev,nconv = eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    v = real.(v)
    println("Eigenvalues: $fs [Hz]")
    
    
    println("Total time elapsed = ",time() - t0,"s")
    
    File =  "fahy_H8.vtk"
    en = 5;
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.VTK.H8;
    scalars=[("Pressure_mode_$en", v[:,en])])
    println("Done")
    true
    
end # fahy_H8_example


function fahy_L2_example()
    println("""
    Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
    Frank J. Fahy, Paolo Gardonio, page 483.
    
    1D mesh.
    """)
    
    t0  =  time()
    
    rho = 1.21*1e-9;# mass density
    c  = 343.0*1000;# millimeters per second
    bulk =  c^2*rho;
    L = 500.0;# length of the box, millimeters
    A = 200.0; # cross-sectional area of the box
    graphics =  true;# plot the solution as it is computed?
    n = 40;#
    neigvs = 8;
    OmegaShift = 10.0;
    
    fens,fes  =  L2block(L,n); # Mesh
    
    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz,1),1))
    
    numberdofs!(P)
    
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(1, 2)), MatAcoustFluid(bulk, rho))
    
    S  =  acousticstiffness(femm, geom, P);
    C  =  acousticmass(femm, geom, P);
    
    d,v,nev,nconv  = eigs(C+OmegaShift*S, S; nev = neigvs, which = :SM)
    d  =  d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    v = real.(v)
    println("Eigenvalues: $fs [Hz]")
    
    
    println("Total time elapsed  =  ",time() - t0,"s")
    
    ix = sortperm(geom.values[:])

    p = PlotlyLight.Plot()
    p(x = geom.values[:][ix], y = v[:,2][ix], type="scatter", mode="lines+markers")
    p.layout.title.text = "Pressure mode 2"
    p.layout.xaxis.title = "x"
    p.layout.yaxis.title = "Pressure Amplitude"
    # plt = lineplot(geom.values[:][ix], v[:,2][ix], title = "Pressure mode 2", canvas = DotCanvas, name = "P", xlabel = "x", ylabel = "P")
    display(p)
    true
    
end # fahy_L2_example

function allrun()
    println("#####################################################") 
    println("# fahy_H20_example ")
    fahy_H20_example()
    println("#####################################################") 
    println("# fahy_H27_example ")
    fahy_H27_example()
    println("#####################################################") 
    println("# fahy_H8_example ")
    fahy_H8_example()
    println("#####################################################") 
    println("# fahy_L2_example ")
    fahy_L2_example()
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module Fahy_examples
