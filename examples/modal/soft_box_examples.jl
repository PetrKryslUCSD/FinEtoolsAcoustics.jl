module soft_box_examples
using FinEtools
using FinEtoolsAcoustics
import Arpack: eigs

function soft_box_Q4_example()
    println("""
    
    Example from Boundary element acoustics: Fundamentals and computer codes, TW Wu, page 123.
    Internal resonance problem. Reference frequencies: 90.7895, 181.579, 215.625, 233.959, 272.368, 281.895
    Quadrilateral mesh.
    """)
    
    t0 = time()
    
    rho=1.21*1e-9;# mass density
    c =345.0*1000;# millimeters per second
    bulk= c^2*rho;
    Lx=1900.0;# length of the box, millimeters
    Ly=800.0; # length of the box, millimeters
    n=14;#
    neigvs=18;
    
    fens,fes = Q4block(Lx,Ly,n,n); # Mesh
    
    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz,1),1))
    bfes = meshboundary(fes)
    setebc!(P, connectednodes(bfes))
    numberdofs!(P)
    
    femm = FEMMAcoust(IntegDomain(fes, GaussRule(2, 2)),
    MatAcoustFluid(bulk,rho))
    
    S = acousticstiffness(femm, geom, P);
    C = acousticmass(femm, geom, P);
    
    d,v,nev,nconv =eigs(C, S; nev=neigvs, which=:SM, explicittransform=:none)
    v = real.(v)
    fs=real(sqrt.(complex(d)))./(2*pi)
    println("Eigenvalues: $fs [Hz]")
    
    
    println("Total time elapsed = ",time() - t0,"s")
    
    File =  "soft_box.vtk"
    scalarllist = Any[]
    for n   in  1:15
        scattersysvec!(P, v[:, n])
        push!(scalarllist, ("Pressure_mode_$n", deepcopy(P.values)));
    end
    vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.VTK.Q4; scalars=scalarllist)
    @async run(`"paraview.exe" $File`)
    
    true
    
end # soft_box_Q4_example

function allrun()
    println("#####################################################") 
    println("# soft_box_Q4_example ")
    soft_box_Q4_example()
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module soft_box_examples
nothing
