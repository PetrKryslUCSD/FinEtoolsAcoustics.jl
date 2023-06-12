
module mmacousticcouplingpanelsm1
using FinEtools
using FinEtoolsAcoustics
using FinEtools.SurfaceNormalModule: SurfaceNormal, updatenormal!
using LinearAlgebra
using Test

function __computenormal!(normalout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt, time = 0.0)
    fill!(normalout, 0.0)
    # We are assuming a surface element here!
    if (size(tangents,1) == 3) && (size(tangents,2) == 2)# surface in three dimensions
        normalout[:] .= cross(vec(tangents[:,1]), vec(tangents[:,2]));# outer normal to the surface
    else
        error("No definition of normal vector");
    end
    nn = norm(normalout);
    if  nn != 0.0 # otherwise return an unnormalized normal
        normalout ./= nn
    end
    return normalout
end

function test()
    R = 0.5 * phun("m")
    E = 205000 * phun("MPa")
    nu = 0.3
    rho = 7850 * phun("kg/m^3");# mass density
    dummybulk, dummyrho = (1.0, 1.0)
    neigvs=6+18;
    OmegaShift=(2*pi*100.0)^2;
    nperradius = 2
    tolerance = R/10000; # geometrical tolerance

    # Construct the mesh
    r(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]]
    origin = [0.0, 0.0, 0.0]
    fens,fes = H8spheren(R, nperradius); # Mesh
    fens1,fes1  =  mirrormesh(fens, fes, [-1.0, 0.0, 0.0], origin,  renumb = r);
    fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1,fes2)
    fens1,fes1  =  mirrormesh(fens, fes, [0.0, -1.0, 0.0], origin,  renumb = r);
    fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1,fes2)
    fens1,fes1  =  mirrormesh(fens, fes, [0.0, 0.0, -1.0], origin,  renumb = r);
    fens,newfes1,fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1,fes2)

    # Extract the mesh of the boundary
    bfes  =  meshboundary(fes)

    # Debugging graphics
    # File  =   "Sphere.vtk"
    # vtkexportmesh(File, fes.conn, fens.xyz, MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    # MR = DeforModelRed3D
    # material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    # femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
    # associategeometry!(femm,  geom) 
    # K  = stiffness(femm, geom, u)
    # M = mass(femm, geom, u)

    # d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    # d = d .- OmegaShift;
    # fs = real(sqrt.(complex(d)))/(2*pi)
    # println("Eigenvalues: $fs [Hz]")

    femm = FEMMAcoustSurf(IntegDomain(bfes, TrapezoidalRule(2)), MatAcoustFluid(dummybulk, dummyrho))
    G = acousticcouplingpanels(femm, geom, u);
    @test abs(G.ff[25 ,   1]  -  0.03125) / abs(0.03125) < 10e-4
    @test abs(G.ff[31 ,   3]  -  0.0255155) / abs(0.0255155) < 10e-4
    @test abs(G.ff[58 ,  13]  -  -0.03125) / abs(-0.03125) < 10e-4
    @test abs(G.ff[91 ,  19]  -  -0.00810978) / abs(-0.00810978) < 10e-4
    @test abs(G.ff[118,  22]  -  0.0105689) / abs(0.0105689) < 10e-4
    @test abs(G.ff[159,  24]  -  -0.019781) / abs(-0.019781) < 10e-4

    G = acousticcouplingpanels(femm, geom, u, SurfaceNormal(3));
    @test abs(G.ff[25 ,   1]  -  0.03125) / abs(0.03125) < 10e-4
    @test abs(G.ff[31 ,   3]  -  0.0255155) / abs(0.0255155) < 10e-4
    @test abs(G.ff[58 ,  13]  -  -0.03125) / abs(-0.03125) < 10e-4
    @test abs(G.ff[91 ,  19]  -  -0.00810978) / abs(-0.00810978) < 10e-4
    @test abs(G.ff[118,  22]  -  0.0105689) / abs(0.0105689) < 10e-4
    @test abs(G.ff[159,  24]  -  -0.019781) / abs(-0.019781) < 10e-4

    G = acousticcouplingpanels(femm, geom, u, SurfaceNormal(3, __computenormal!));
    @test abs(G.ff[25 ,   1]  -  0.03125) / abs(0.03125) < 10e-4
    @test abs(G.ff[31 ,   3]  -  0.0255155) / abs(0.0255155) < 10e-4
    @test abs(G.ff[58 ,  13]  -  -0.03125) / abs(-0.03125) < 10e-4
    @test abs(G.ff[91 ,  19]  -  -0.00810978) / abs(-0.00810978) < 10e-4
    @test abs(G.ff[118,  22]  -  0.0105689) / abs(0.0105689) < 10e-4
    @test abs(G.ff[159,  24]  -  -0.019781) / abs(-0.019781) < 10e-4



    true

end
end
using .mmacousticcouplingpanelsm1
mmacousticcouplingpanelsm1.test()
