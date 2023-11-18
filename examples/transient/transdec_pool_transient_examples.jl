module transdec_pool_transient_examples
using FinEtools
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection, vtkwrite
using FinEtoolsAcoustics
using LinearAlgebra


function transdec_pool_transient(freq = 100)
	rho = 1000.0*phun("kg/m^3");# mass density
        c  = 1500.0*phun("m/s");# sound speed
	bulk =  c^2*rho;
	omega =  freq*phun("rev/s");      # frequency of the piston
	P_piston = 1.0e-3*phun("MPa"); # amplitude of the piston pressure
	R = 0.5*phun("ft");# radius of the piston
    rim_distance = 150 * phun("ft")
	tolerance = R/100
	dt = 1.0/freq/20;
	tfinal = rim_distance / c  * 3 / 5
	nsteps = Int(round(tfinal/dt))+1;
    nbtw = max(1, Int(round(nsteps/300)))
    pos_pressure_fun(t) = P_piston*sin(omega*t)
    pos_pressuredd_fun(t) = (-omega^2)*pos_pressure_fun(t)
    neg_pressure_fun(t) = P_piston*sin(omega*t)
    neg_pressuredd_fun(t) = (-omega^2)*neg_pressure_fun(t)

    t0  =  time()

    # Hexahedral mesh
    mesh = import_ABAQUS(joinpath(@__DIR__, "transdec-pool-quarter.inp"))
    fens = mesh["fens"]
    fes = mesh["fesets"][1]

    fens.xyz *= phun("ft")

    @info "$(count(fens)) nodes"

    fens1, fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 25.0, 56.65].*phun("ft"); renumb =  c ->  c[[1, 3, 2, 4]])
    fens, newfes1, fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    bfes  =  meshboundary(fes)
    File  =   "transdec-pool-boundary.vtk"
    vtkexportmesh(File, fens, bfes)
    # @async run(`"paraview.exe" $File`)

    l1 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])
    l2 = selectelem(fens, bfes, facing = true, direction = [+1.0 0.0 0.0])
    l3 = selectelem(fens, bfes, distance = R, from = [0.0, 6.0, 56.65].*phun("ft"), inflate = tolerance)
    l4 = selectelem(fens, bfes, box = [-Inf, Inf, 25.0, 25.0, -Inf, Inf].*phun("ft"), inflate = tolerance)

    lpos = intersect(l2, l3)
    lneg = intersect(l1, l3)

    File  =   "transdec-pool-pos.vtk"
    vtkexportmesh(File, fens, subset(bfes, lpos))

    File  =   "transdec-pool-neg.vtk"
    vtkexportmesh(File, fens, subset(bfes, lneg))

    File  =   "transdec-pool-free.vtk"
    vtkexportmesh(File, fens, subset(bfes, l4))

    # Piston surface mesh
    piston_pos_fes = subset(bfes, lpos);
    piston_neg_fes = subset(bfes, lneg);

    println("Pre-processing time elapsed  =  ",time() - t0,"s")

	t1  =  time()

	material = MatAcoustFluid(bulk, rho)
	# Region of the fluid
	femm = FEMMAcoust(IntegDomain(fes, TetRule(1)), material)

	geom = NodalField(fens.xyz);
	P = NodalField(zeros(nnodes(geom),1));
	piston_pos_fenids = connectednodes(piston_pos_fes);
	setebc!(P, piston_pos_fenids, true, 1, 0.0);
    piston_neg_fenids = connectednodes(piston_neg_fes);
    setebc!(P, piston_neg_fenids, true, 1, 0.0);
    setebc!(P, connectednodes(subset(bfes, l4)), true, 1, 0.0);
	applyebc!(P);
	numberdofs!(P);

	S  =  acousticstiffness(femm, geom, P);
	C  =  acousticmass(femm, geom, P);
	D= (2.0/dt)*S +(dt/2.0)*C;

	P0 = deepcopy(P)
	Pdd0 = deepcopy(P)
	P1 = deepcopy(P)
	Pdd1 = deepcopy(P)
	TMPF = deepcopy(P)
	P = nothing # we don't need this field anymore
	vP0 = zeros(P0.nfreedofs)
	vP1 = deepcopy(vP0)
	vPd0 = deepcopy(vP0)
	vPd1 = deepcopy(vP0)

	t = 0.0
	P0.fixed_values[piston_pos_fenids,1] .= pos_pressure_fun(t)
	Pdd0.fixed_values[piston_pos_fenids,1] .= pos_pressuredd_fun(t)
    P0.fixed_values[piston_neg_fenids,1] .= neg_pressure_fun(t)
    Pdd0.fixed_values[piston_neg_fenids,1] .= neg_pressuredd_fun(t)
	vP0 = gathersysvec!(P0, vP0)

	# nh = selectnode(fens, nearestto = [R+Ro/2, 0.0, 0.0] )
	# Pnh = [P1.values[nh, 1][1]]
	pressures = [(t, vP1)]
	step = 0;
	while t <=tfinal
		step = step  +1;
		t=t+dt;
		P1.fixed_values[piston_pos_fenids,1] .= pos_pressure_fun(t)
		Pdd1.fixed_values[piston_pos_fenids,1] .= pos_pressuredd_fun(t)
        P1.fixed_values[piston_neg_fenids,1] .= neg_pressure_fun(t)
        Pdd1.fixed_values[piston_neg_fenids,1] .= neg_pressuredd_fun(t)
		TMPF.fixed_values = P0.fixed_values + P1.fixed_values
		F = nzebcloadsacousticmass(femm, geom, TMPF);
		TMPF.fixed_values = Pdd0.fixed_values + Pdd1.fixed_values
		F = F + nzebcloadsacousticstiffness(femm, geom, TMPF);
		# println("$(norm(F))")
		vPd1 = D\((2/dt)*(S*vPd0) - C*(2*vP0+(dt/2)*vPd0) + F);
		vP1 = vP0 + (dt/2)*(vPd0+vPd1);
		scattersysvec!(P1, vP1); # store current pressure
        if rem(step+1, nbtw) == 0
            push!(pressures, (t, deepcopy(vP1)))
        end
        # Swap variables for the next step
		copyto!(vP0, vP1)
		copyto!(vPd0, vPd1)
		copyto!(P0, P1)
		copyto!(Pdd0, Pdd1)
		println("step = $( step )/$(nsteps)")
	end

    # Visualization
    @info "Dumping visualization"
    times = Float64[]
    scalars = []
    for i in 1:length(pressures)
        t, p = pressures[i]
        scattersysvec!(P1, p)
        push!(scalars, ("p", deepcopy(P1.values[:, 1])))
        push!(times, t)
    end
    vtkwritecollection("transdec_pool_transient-p", fens, fes, times; scalars = scalars)

	true
end #

function allrun()
    println("#####################################################") 
    println("# transdec_pool_example ")
    transdec_pool_transient(100)
    # transdec_pool_transient(200)
    # transdec_pool_transient(500)
    # transdec_pool_transient(750)
    # transdec_pool_transient(1000)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module baffled_piston_examples
nothing
