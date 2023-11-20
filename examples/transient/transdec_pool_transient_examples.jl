module transdec_pool_transient_examples
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtools.MeshExportModule.VTKWrite: vtkwritecollection, vtkwrite
using FinEtoolsAcoustics
using DataDrop
using LinearAlgebra
using PlotlyLight

function transdec_pool_harmonic(freq = 100)
	rho = 1000.0*phun("kg/m^3");# mass density
        c  = 1500.0*phun("m/s");# sound speed
	bulk =  c^2*rho;
	omega =  freq*phun("rev/s");      # frequency of the piston
	P_piston = 1.0e-3*phun("MPa"); # amplitude of the piston pressure
	R = 0.5*phun("ft");# radius of the piston
    rim_distance = 150 * phun("ft")
	tolerance = R/100
	dt = 1.0/freq/20;
	tfinal = rim_distance / c  * 3
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
	C = acousticmass(femm, geom, P);
	D = (2.0/dt)*S +(dt/2.0)*C;

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
	P0.values[piston_pos_fenids,1] .= pos_pressure_fun(t)
	Pdd0.values[piston_pos_fenids,1] .= pos_pressuredd_fun(t)
    P0.values[piston_neg_fenids,1] .= neg_pressure_fun(t)
    Pdd0.values[piston_neg_fenids,1] .= neg_pressuredd_fun(t)
	vP0 = gathersysvec!(P0, vP0)

	# nh = selectnode(fens, nearestto = [R+Ro/2, 0.0, 0.0] )
	# Pnh = [P1.values[nh, 1][1]]
	pressures = [(t, vP1)]
	step = 0;
	while t <=tfinal
		step = step  +1;
		t=t+dt;
		P1.values[piston_pos_fenids,1] .= pos_pressure_fun(t)
		Pdd1.values[piston_pos_fenids,1] .= pos_pressuredd_fun(t)
        P1.values[piston_neg_fenids,1] .= neg_pressure_fun(t)
        Pdd1.values[piston_neg_fenids,1] .= neg_pressuredd_fun(t)
		TMPF.values = P0.values + P1.values
		F = nzebcloadsacousticmass(femm, geom, TMPF);
		TMPF.values = Pdd0.values + Pdd1.values
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
    vtkwritecollection("transdec_pool_harmonic-f=$(freq)-p", fens, fes, times; scalars = scalars)

	true
end #

function transect_interpolation(fens, fes, from, to, npoints, geometricaltolerance)
    parametrictolerance = 1.0e-5
    distances = Float64[]
    nodebox = initbox!([], vec(from))
    for i in 1:npoints
        p = (i - npoints) / (1 - npoints) .* from + (i - 1) / (npoints - 1) .* to
        push!(distances, norm(p - from))
        updatebox!(nodebox, vec(p))
    end
    nodebox = inflatebox!(nodebox, geometricaltolerance)
    el = selectelem(fens, fes; overlappingbox = nodebox)
    interp = []
    for i in 1:npoints
        p = (i - npoints) / (1 - npoints) .* from + (i - 1) / (npoints - 1) .* to
        for e in el
            c = [k for k in fes.conn[e]]
            pc, success = map2parametric(fes,
                fens.xyz[c, :],
                vec(p);
                tolerance = parametrictolerance,
                maxiter = 7,)
            @assert success # this shouldn't be tripped; normally we succeed
            if inparametric(fes, pc; tolerance = parametrictolerance)
                N = bfun(fes, pc)
                push!(interp, (c, N))
                # ff.values[i, :] = transpose(N) * fcsub.values[c, :]
                break
            end
        end
    end
    @assert length(interp) == npoints
    return distances, interp
end

function transdec_pool_pulse(freq = 100)
    rho = 1000.0*phun("kg/m^3");# mass density
        c  = 1500.0*phun("m/s");# sound speed
    bulk =  c^2*rho;
    omega =  freq*phun("rev/s");      # frequency of the piston
    P_piston = 1.0e-3*phun("MPa"); # amplitude of the piston pressure
    R = 0.5*phun("ft");# radius of the piston
    rim_distance = 150 * phun("ft")
    tolerance = R/100
    dt = 1.0/freq/20;
    tfinal = rim_distance / c  * 0.2
    nsteps = Int(round(tfinal/dt))+1;
    nbtw = max(1, Int(round(nsteps/300)))
    pos_pressure_fun(t) = (t <= 1/freq) ? P_piston*sin(omega*t) : 0.0
    pos_pressuredd_fun(t) = (-omega^2)*pos_pressure_fun(t)
    neg_pressure_fun(t) = (t <= 1/freq) ? P_piston*sin(omega*t) : 0.0
    neg_pressuredd_fun(t) = (-omega^2)*neg_pressure_fun(t)

    t0 = time()

    # Hexahedral mesh
    mesh = import_ABAQUS(joinpath(@__DIR__, "transdec-pool-quarter.inp"))
    fens = mesh["fens"]
    fes = mesh["fesets"][1]

    fens.xyz *= phun("ft")

    @info "$(count(fens)) nodes"

    fens1, fes1 = mirrormesh(fens, fes, [-1.0, 0.0, 0.0], [0.0, 25.0, 56.65].*phun("ft"); renumb =  c ->  c[[1, 3, 2, 4]])
    fens, newfes1, fes2 =  mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    bfes = meshboundary(fes)
    File =  "transdec-pool-boundary.vtk"
    vtkexportmesh(File, fens, bfes)
    # @async run(`"paraview.exe" $File`)

    l1 = selectelem(fens, bfes, facing = true, direction = [-1.0 0.0 0.0])
    l2 = selectelem(fens, bfes, facing = true, direction = [+1.0 0.0 0.0])
    l3 = selectelem(fens, bfes, distance = R, from = [0.0, 6.0, 56.65].*phun("ft"), inflate = tolerance)
    l4 = selectelem(fens, bfes, box = [-Inf, Inf, 25.0, 25.0, -Inf, Inf].*phun("ft"), inflate = tolerance)

    lpos = intersect(l2, l3)
    lneg = intersect(l1, l3)

    File =  "transdec-pool-pos.vtk"
    vtkexportmesh(File, fens, subset(bfes, lpos))

    File =  "transdec-pool-neg.vtk"
    vtkexportmesh(File, fens, subset(bfes, lneg))

    File =  "transdec-pool-free.vtk"
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

    S_ff, S_fd = matrix_blocked(S, nfreedofs(P))
    C_ff, C_fd = matrix_blocked(C, nfreedofs(P))
    D_ff = (2.0/dt)*S_ff +(dt/2.0)*C_ff;

    P0 = deepcopy(P)
    Pdd0 = deepcopy(P)
    P1 = deepcopy(P)
    Pdd1 = deepcopy(P)
    TMPF = deepcopy(P)
    P = nothing # we don't need this field anymore
    vP0 = gathersysvec(P0, :f)
    vP1 = deepcopy(vP0)
    vPd0 = deepcopy(vP0)
    vPd1 = deepcopy(vP0)
    F_f = deepcopy(vP0)
    vTMP_d = gathersysvec(P0, :d)

    t = 0.0
    P0.values[piston_pos_fenids,1] .= pos_pressure_fun(t)
    Pdd0.values[piston_pos_fenids,1] .= pos_pressuredd_fun(t)
    P0.values[piston_neg_fenids,1] .= neg_pressure_fun(t)
    Pdd0.values[piston_neg_fenids,1] .= neg_pressuredd_fun(t)
    vP0 = gathersysvec!(P0, vP0)

    # nh = selectnode(fens, nearestto = [R+Ro/2, 0.0, 0.0] )
    # Pnh = [P1.values[nh, 1][1]]
    pressures = [(t, deepcopy(P1.values))]
    step = 0;
    while t <=tfinal
        step = step  +1;
        t=t+dt;
        P1.values[piston_pos_fenids,1] .= pos_pressure_fun(t)
        Pdd1.values[piston_pos_fenids,1] .= pos_pressuredd_fun(t)
        P1.values[piston_neg_fenids,1] .= neg_pressure_fun(t)
        Pdd1.values[piston_neg_fenids,1] .= neg_pressuredd_fun(t)
        TMPF.values = P0.values + P1.values
        F_f .= - C_fd * gathersysvec!(TMPF, vTMP_d, :d)
        TMPF.values = Pdd0.values + Pdd1.values
        F_f .+= - S_fd * gathersysvec!(TMPF, vTMP_d, :d)
        vPd1 = D_ff \ ((2/dt)*(S_ff*vPd0) - C_ff*(2*vP0+(dt/2)*vPd0) + F_f);
        vP1 = vP0 + (dt/2)*(vPd0+vPd1);
        scattersysvec!(P1, vP1); # store current pressure
        if rem(step+1, nbtw) == 0
            push!(pressures, (t, deepcopy(P1.values)))
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

    from = [-0.01524000021336, 1.8288000000000002, 17.2669204572]
    to = from .+ [-1.5, 0, 0]
    npoints = 10
    geometricaltolerance = tolerance / 100
    distances, interp = transect_interpolation(fens, fes, from, to, npoints, geometricaltolerance)

    # plot = PlotlyLight.Plot()
    # plot.layout.yaxis.title = "Distance"
    # plot.layout.xaxis.title = "Pressure"
    # plot.layout.xaxis.side = "top"
    # plot.layout.margin.pad = 10

    times = Float64[]
    scalars = []
    transects = []
    for i in 1:length(pressures)
        t, p = pressures[i]
        push!(scalars, ("p", deepcopy(p)))
        ptrans = zeros(npoints)
        for (j, ip) in enumerate(interp)
            ptrans[j] = dot(ip[2], p[ip[1]])
        end
        push!(transects, ptrans)
        # plot(x = vec(distances), y = vec(ptrans), mode="lines")
        # display(plot)
        push!(times, t)
    end
    vtkwritecollection("transdec_pool_pulse-f=$(freq)-p", fens, fes, times; scalars = scalars)
    d = Dict("times" => times,
            "distances" => distances,
            "transects" => transects
        )
    DataDrop.store_json("transdec_pool_pulse-f=$(freq)-ptrans", d)

    true
end #

function allrun()
    println("#####################################################") 
    println("# transdec_pool_example ")
    # transdec_pool_harmonic(100)
    # transdec_pool_harmonic(200)
    # transdec_pool_harmonic(500)
    # transdec_pool_harmonic(750)
    # transdec_pool_harmonic(1000)
    transdec_pool_pulse(100)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module transdec_pool_transient_examples

nothing
