module sphere_dipole_examples
using FinEtools
using FinEtoolsAcoustics
using LinearAlgebra


function sphere_dipole_transient()
    println("The interior sphere accelerates in the alternately in the positive
    and negative x-direction, generating positive pressure ahead of it, negative
    pressure behind. Time-dependent simulation.
    ")
    rho = 1.21*phun("kg/m^3");# mass density
    c  = 343.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    a_amplitude=1.0*phun("mm/s^2");# amplitude of the  acceleration of the sphere
    # omega = 2000*phun("rev/s");      # frequency of the incident wave
    R = 50.0*phun("mm"); # radius of the interior sphere
    Ro = 8*R # radius of the external sphere
    P_amplitude = R*rho*a_amplitude; # pressure amplitude
    nref=2;
    nlayers=40;
    tolerance = Ro/(nlayers)/100
    frequency=1200.; # Hz
    omega=2*pi*frequency;
    dt=1.0/frequency/20;
    tfinal=90*dt;
    nsteps=round(tfinal/dt)+1;
    
    t0  =  time()
    
    # Hexahedral mesh
    fens,fes  =  H8sphere(R,nref);
    fens1,fes1  =  mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0],
    renumb =  r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]]);
    fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1,fes2)
    
    
    # Derive the finite element sets for the boundary
    bfes  =  meshboundary(fes)
    # Outer spherical boundary
    function dout(xyz)
        return xyz/norm(xyz)
    end
    
    louter = selectelem(fens, bfes, facing = true, direction = dout)
    outer_fes = subset(bfes, louter);
    
    fens,fes  = H8extrudeQ4(fens, outer_fes, nlayers, (xyz, layer)->(R+layer/nlayers*(Ro-R))*xyz/norm(xyz));
    connected = findunconnnodes(fens, fes); fens, new_numbering = compactnodes(fens, connected);
    fess = renumberconn!(fes, new_numbering);
    
    geom  =  NodalField(fens.xyz)
    P  =  NodalField(zeros(FCplxFlt,size(fens.xyz,1),1))
    
    bfes  =  meshboundary(fes)
    
    # File  =   "Sphere.vtk"
    # vtkexportmesh(File, bfes.conn, geom.values, FinEtools.MeshExportModule.Q4)
    # @async run(`"paraview.exe" $File`)
    
    # File  =   "Sphere.vtk"
    # vtkexportmesh (File, fes.conn, fens.xyz, MeshExportModule.H8)
    # @async run(`"C:/Program Files (x86)/ParaView 4.2.0/bin/paraview.exe" $File`)
    
    linner = selectelem(fens, bfes, distance = R, from = [0.0 0.0 0.0], inflate = tolerance)
    louter = selectelem(fens, bfes, facing = true, direction = dout)
    
    println("Pre-processing time elapsed  =  ",time() - t0,"s")
    
    t1  =  time()
    
    numberdofs!(P)
    
    material = MatAcoustFluid(bulk,rho)
    femm  =  FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), material)
    
    S  =  acousticstiffness(femm, geom, P);
    C  =  acousticmass(femm, geom, P);
    
    abcfemm  =  FEMMAcoustSurf(IntegDomain(subset(bfes, louter), GaussRule(2, 2)), material)
    D  =  acousticABC(abcfemm, geom, P);
    
    # Inner sphere pressure loading
    function dipole(dpdn, xyz, J, label, t)
        n = cross(J[:,1],J[:,2]);
        n = vec(n/norm(n));
        dpdn[1] = -rho*a_amplitude*sin(omega*t)*n[1]
    end
    
    dipfemm  =  FEMMAcoustSurf(IntegDomain(subset(bfes, linner), GaussRule(2, 2)), material)

    # Solve
    P0 = deepcopy(P)
    P0.values[:] .= 0.0; # initially all pressure is zero
    vP0 = gathersysvec(P0);
    vP1 = zeros(eltype(vP0), size(vP0));
    vQ0 = zeros(eltype(vP0), size(vP0));
    vQ1 = zeros(eltype(vP0), size(vP0));
    t = 0.0;
    P1 = deepcopy(P0);
    
    fi  =  ForceIntensity(FCplxFlt, 1, (dpdn, xyz, J, label)->dipole(dpdn, xyz, J, label, t));
    La0 = distribloads(dipfemm, geom, P1, fi, 2);
    
    A = (2.0/dt)*S + D + (dt/2.)*C;
    
    step =0;
    while t <= tfinal
      step = step  + 1;
      println("Time $t ($(step)/$(round(tfinal/dt)+1))")
      t = t+dt;
      fi  = ForceIntensity(FCplxFlt, 1, (dpdn, xyz, J, label)->dipole(dpdn, xyz, J, label, t));
      La1 = distribloads(dipfemm, geom, P1, fi, 2);
      vQ1 = A\((2/dt)*(S*vQ0)-D*vQ0-C*(2*vP0+(dt/2)*vQ0)+La0+La1);
      vP1 = vP0 + (dt/2)*(vQ0+vQ1);
      vP0 = vP1;
      vQ0 = vQ1;
      P1 = scattersysvec!(P1, vec(vP1));
      P0 = P1;
      La0 = La1;
    end
    
    File  =   "sphere_dipole_1.vtk"
    vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.VTK.H8;
     scalars = [( "realP", real(P1.values))])
    @async run(`"paraview.exe" $File`)
    
    
    true
end # sphere_dipole_transient

function allrun()
    println("#####################################################") 
    println("# sphere_dipole_transient ")
    sphere_dipole_transient()
    return true
end # function allrun

end # module baffled_piston_examples
