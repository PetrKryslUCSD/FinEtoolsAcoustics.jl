module straight_duct_examples
using FinEtools
using FinEtoolsAcoustics
using PlotlyLight

function straight_duct_Q8_example()
    t0  =  time()
    
    rho = 1.21*phun("kg/m^3");# mass density
    c  = 343.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    omega =  54.5901*phun("rev/s")
    vn0 =  -1.0*phun("m/s")
    Lx = 10.0*phun("m");# length of the box, millimeters
    Ly = 1.0*phun("m"); # length of the box, millimeters
    n = 20;#number of elements along the length
    
    println("""
    
    Straight duct with anechoic termination.
    Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 44.
    Both real and imaginary components of the pressure should have amplitude of
    rho*c = $(rho*c).
    
    Two-dimensional model. Quadratic qquadrilateral mesh.
    """)
    
    fens,fes  =  Q8block(Lx,Ly,n,2); # Mesh
    bfes  =  meshboundary(fes)
    L0 = selectelem(fens,bfes,facing = true, direction = [-1.0 0.0])
    L10 = selectelem(fens,bfes,facing = true, direction = [+1.0 0.0])
    nLx = selectnode(fens,box = [0.0 Lx  0.0 0.0], inflate = Lx/1.0e5)
    
    geom  =  NodalField(fens.xyz)
    P  =  NodalField(zeros(ComplexF64,size(fens.xyz,1),1))
    
    numberdofs!(P)
    
    
    material = MatAcoustFluid(bulk,rho)
    femm  =  FEMMAcoust(IntegDomain(fes, GaussRule(2, 3)), material)
    
    S  =  shouldbemass(femm, geom, P);
    C  =  shouldbestiffness(femm, geom, P);
    
    
    E10femm  =  FEMMAcoustSurf(IntegDomain(subset(bfes,L10),GaussRule(1, 3)), material)
    D  =  acousticABC(E10femm, geom, P);
    
    E0femm  =  FEMMBase(IntegDomain(subset(bfes,L0), GaussRule(1,  3)))
    fi  =  ForceIntensity(-1.0im*omega*rho*vn0);
    F  =  distribloads(E0femm, geom, P, fi, 2);
    
    p = (-omega^2*S +omega*1.0im*D + C)\F
    scattersysvec!(P, p[:])
    
    println("Pressure amplitude bounds: ")
    println("  real $(minimum(real(P.values)))/$(maximum(real(P.values)))")
    println("  imag $(minimum(imag(P.values)))/$(maximum(imag(P.values)))")
    
    println("Total time elapsed  =  ",time() - t0,"s")
    
    File  =   "straight_duct_Q8.vtk"
    scalars = real(P.values);
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.VTK.Q8;  scalars = [("Pressure", scalars)])
    @async run(`"paraview.exe" $File`)
    
    ix = sortperm(geom.values[nLx,1])

    p = PlotlyLight.Plot()
    p(x = geom.values[nLx,1][ix], y = real(P.values)[nLx][ix], type="scatter", mode="lines+markers")
    p.layout.title.text = "Steady-state pressure"
    p.layout.xaxis.title = "x"
    p.layout.yaxis.title = "Pressure Amplitude"

    display(p)

    true
    
end # straight_duct_Q8_example


function straight_duct_T3_example()
    t0  =  time()
    
    rho = 1.21*phun("kg/m^3");# mass density
    c  = 343.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    omega =  54.5901*phun("rev/s")
    vn0 =  -1.0*phun("m/s")
    Lx = 10.0*phun("m");# length of the box, millimeters
    Ly = 1.0*phun("m"); # length of the box, millimeters
    n = 20;#number of elements along the length
    
    println("""
    
    Straight duct with anechoic termination.
    Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 44.
    Both real and imaginary components of the pressure should have amplitude of
    rho*c = $(rho*c).
    
    Triangle mesh.
    """)
    
    fens,fes  =  T3block(Lx,Ly,n,2); # Mesh
    bfes  =  meshboundary(fes)
    L0 = selectelem(fens,bfes,facing = true, direction = [-1.0 0.0])
    L10 = selectelem(fens,bfes,facing = true, direction = [+1.0 0.0])
    nLx = selectnode(fens,box = [0.0 Lx  0.0 0.0], inflate = Lx/1.0e5)
    
    geom  =  NodalField(fens.xyz)
    P  =  NodalField(zeros(ComplexF64,size(fens.xyz,1),1))
    
    numberdofs!(P)
    
    material = MatAcoustFluid(bulk,rho)
    femm  =  FEMMAcoust(IntegDomain(fes, TriRule(1)), material)
    
    S  =  shouldbemass(femm, geom, P);
    C  =  shouldbestiffness(femm, geom, P);
    
    
    E10femm  =  FEMMAcoustSurf(IntegDomain(subset(bfes,L10),GaussRule(1, 2)), material)
    D  =  acousticABC(E10femm, geom, P);
    
    E0femm  =  FEMMBase(IntegDomain(subset(bfes,L0), GaussRule(1,  2)))
    fi  =  ForceIntensity(-1.0im*omega*rho*vn0);
    F  =  distribloads(E0femm, geom, P, fi, 2);
    
    p = (-omega^2*S +omega*1.0im*D + C)\F
    scattersysvec!(P, p[:])
    
    println("Pressure amplitude bounds: ")
    println("  real $(minimum(real(P.values)))/$(maximum(real(P.values)))")
    println("  imag $(minimum(imag(P.values)))/$(maximum(imag(P.values)))")
    
    println("Total time elapsed  =  ",time() - t0,"s")
    
    File  =   "straight_duct.vtk"
    scalars = real(P.values);
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.VTK.T3; scalars = [("Pressure", scalars)])
    @async run(`"paraview.exe" $File`)
     
     ix = sortperm(geom.values[nLx,1])

     p = PlotlyLight.Plot()
     p(x = geom.values[nLx,1][ix], y = real(P.values)[nLx][ix], type="scatter", mode="lines+markers")
     p.layout.title.text = "Steady-state pressure"
     p.layout.xaxis.title = "x"
     p.layout.yaxis.title = "Pressure Amplitude"

     display(p)
    true
    
end # straight_duct_T3_example

function allrun()
    println("#####################################################") 
    println("# straight_duct_Q8_example ")
    straight_duct_Q8_example()
    println("#####################################################") 
    println("# straight_duct_T3_example ")
    straight_duct_T3_example()
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module straight_duct_examples
nothing
