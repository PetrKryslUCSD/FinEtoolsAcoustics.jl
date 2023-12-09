
module mmfahyH8example2
using LinearAlgebra
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
P = NodalField(fill(zero(FFlt), size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoust(IntegDomain(fes, GaussRule(3, 2)), MatAcoustFluid(bulk, rho))


Ma = acousticmass(femm, geom, P);
Ka = acousticstiffness(femm, geom, P);

d,v,nev,nconv = eigs(Ka+OmegaShift*Ma, Ma; nev=neigvs, which=:SM)
d = d .- OmegaShift;
fs = real(sqrt.(complex(d)))/(2*pi)
# println("Eigenvalues: $fs [Hz]")
#
#
# println("Total time elapsed = ",time() - t0,"s")
@test norm([0.00000e+00, 3.43088e+02, 6.86705e+02, 1.03138e+03, 1.37765e+03, 1.72604e+03, 
2.07709e+03, 2.43134e+03]  .- fs) ./ norm(fs)   < 1.0e-3
# File =  "fahy_H8.vtk"
# en = 5;
# vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
#           scalars=[("Pressure_mode_$en", v[:,en])])
# println("Done")

  @test (fs[2]-343.08)/343.08 < 0.001

  true
end
end
using .mmfahyH8example2
mmfahyH8example2.test()

module mmfahyT4example
using LinearAlgebra
using FinEtools
using FinEtoolsAcoustics
using Test
import Arpack: eigs
function test()

# println("""
# Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
# Frank J. Fahy, Paolo Gardonio, page 483.
#
# Tetrahedral mesh.
# """)

t0 = time()

rho=1.21*1e-9;# mass density
c =343.0*1000;# millimeters per second
bulk= c^2*rho;
L=500.0;# length of the box, millimeters
A=200.0; # cross-sectional area of the box
n=40;#
neigvs=8;
OmegaShift=10.0;

fens,fes = T4block(L,sqrt(A),sqrt(A),n,1,1); # Mesh

geom = NodalField(fens.xyz)
P = NodalField(fill(zero(FFlt), size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoust(IntegDomain(fes, TetRule(1)), MatAcoustFluid(bulk, rho))


Ma = acousticmass(femm, geom, P);
Ka = acousticstiffness(femm, geom, P);

d,v,nev,nconv = eigs(Ka+OmegaShift*Ma, Ma; nev=neigvs, which=:SM)
d = d .- OmegaShift;
fs = real(sqrt.(complex(d)))/(2*pi)
# println("Eigenvalues: $fs [Hz]")

@test norm([0.00000e+00, 3.43131e+02, 6.87046e+02, 1.03253e+03, 1.38034e+03, 1.73126e+03, 2.08602e+03, 2.44538e+03]   .- fs) ./ norm(fs)   < 1.0e-3
#
#
# println("Total time elapsed = ",time() - t0,"s")

# File =  "fahy_H8.vtk"
# en = 5;
# vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
#           scalars=[("Pressure_mode_$en", v[:,en])])
# println("Done")

  @test (fs[2]-343.08)/343.08 < 0.001

  true
end
end
using .mmfahyT4example
mmfahyT4example.test()


module mmfahyT4NICEexample
using LinearAlgebra
using FinEtools
using FinEtoolsAcoustics
using Test
import Arpack: eigs
function test()

# println("""
# Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
# Frank J. Fahy, Paolo Gardonio, page 483.
#
# Tetrahedral NICE mesh.
# """)

t0 = time()

rho=1.21*1e-9;# mass density
c =343.0*1000;# millimeters per second
bulk= c^2*rho;
L=500.0;# length of the box, millimeters
A=200.0; # cross-sectional area of the box
n=40;#
neigvs=8;
OmegaShift=10.0;

fens,fes = T4block(L,sqrt(A),sqrt(A),n,1,1); # Mesh

geom = NodalField(fens.xyz)
P = NodalField(fill(zero(FFlt), size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoust(IntegDomain(fes, TetRule(1)), MatAcoustFluid(bulk, rho))
S = acousticmass(femm, geom, P);
femm = FEMMAcoustNICE(IntegDomain(fes, NodalSimplexRule(3)), MatAcoustFluid(bulk, rho))
associategeometry!(femm,  geom)
C = acousticstiffness(femm, geom, P);

d,v,nev,nconv = eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
d = d .- OmegaShift;
fs = real(sqrt.(complex(d)))/(2*pi)
# println("Eigenvalues: $fs [Hz]")

@test norm([5.13447e-05, 3.42904e+02, 6.85231e+02, 1.02641e+03, 1.36585e+03, 1.70298e+03, 2.03722e+03, 2.36798e+03]  .- fs) ./ norm(fs)   < 1.0e-3
#
#
# println("Total time elapsed = ",time() - t0,"s")

# File =  "fahy_H8.vtk"
# en = 5;
# vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
#           scalars=[("Pressure_mode_$en", v[:,en])])
# println("Done")

  @test (fs[2]-343.08)/343.08 < 0.001

  true
end
end
using .mmfahyT4NICEexample
mmfahyT4NICEexample.test()


module mmfahyT4NICEexample2
using LinearAlgebra
using FinEtools
using FinEtoolsAcoustics
using Test
import Arpack: eigs
function test()

# println("""
# Example from Sound and Structural Vibration, Second Edition: Radiation, Transmission and Response [Paperback]
# Frank J. Fahy, Paolo Gardonio, page 483.
#
# Tetrahedral NICE mesh.
# """)

t0 = time()

rho=1.21*1e-9;# mass density
c =343.0*1000;# millimeters per second
bulk= c^2*rho;
L=500.0;# length of the box, millimeters
A=200.0; # cross-sectional area of the box
n=80;#
neigvs=8;
OmegaShift=10.0;

fens,fes = T4block(L,sqrt(A),sqrt(A), n, 6, 6); # Mesh

geom = NodalField(fens.xyz)
P = NodalField(fill(zero(FFlt), size(fens.xyz,1),1))

numberdofs!(P)

femm = FEMMAcoust(IntegDomain(fes, TetRule(1)), MatAcoustFluid(bulk, rho))
S = acousticmass(femm, geom, P);
femm = FEMMAcoustNICE(IntegDomain(fes, NodalSimplexRule(3)), MatAcoustFluid(bulk, rho))
associategeometry!(femm,  geom)
C = acousticstiffness(femm, geom, P);

d,v,nev,nconv = eigs(C+OmegaShift*S, S; nev=neigvs, which=:SM)
d = d .- OmegaShift;
fs = real(sqrt.(complex(d)))/(2*pi)
# println("Eigenvalues: $fs [Hz]")

@test norm([0.00000e+00, 3.42968e+02, 6.85747e+02, 1.02815e+03, 1.36997e+03, 1.71104e+03, 2.05116e+03, 2.39013e+03]  .- fs) ./ norm(fs)   < 1.0e-3
#
#
# println("Total time elapsed = ",time() - t0,"s")

# File =  "fahy_H8.vtk"
# en = 5;
# vtkexportmesh(File, fes.conn, geom.values, FinEtools.MeshExportModule.H8;
#           scalars=[("Pressure_mode_$en", v[:,en])])
# println("Done")

  @test (fs[2]-343.08)/343.08 < 0.001

  true
end
end
using .mmfahyT4NICEexample2
mmfahyT4NICEexample2.test()
