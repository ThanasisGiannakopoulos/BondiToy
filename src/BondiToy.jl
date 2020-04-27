module BondiToy

export Param, run_toy

using Parameters
using DifferentialEquations
using Interpolations
using HDF5
using DelimitedFiles
using Random
using Printf

@with_kw struct Param
    # discretization parameters
    NX                :: Int
    Nz                :: Int
    cX                :: Float64 = 10.0
    rmin              :: Float64 = 2.0
    umax              :: Float64
    # parameter that controls amplitude of noise
    noise_amplitude   :: Float64
    # directory to save data
    out_dir           :: String
    out_every         :: Int
    compute_L2_norm   :: Bool = true
    compute_Lop_norm  :: Bool = true
end

rtoX(r, cX, rmin) = (r .- rmin) ./ sqrt.(cX*cX .+ (r .- rmin).^2)
Xtor(X, cX, rmin) = rmin .+ cX * X ./ sqrt.(1.0 .- X.^2)
dX_dr(X, cX) = (1.0 .- X.^2).^1.5 ./ cX

struct System
    X    :: Vector
    r    :: Vector
    dXdr :: Vector
    hX   :: Float64
    z    :: Vector
    hz   :: Float64
end
function System(p::Param)
    hX = 1.0 / (p.NX-1)
    X  = range(0, 1, length=p.NX)

    hz = 2.0*π / (p.Nz) # remove last point; 0=2π (periodic)
    z  = range(0, 2.0*π-hz, length=p.Nz)

    rr    = Xtor(X, p.cX, p.rmin)
    dXdr  = dX_dr(X, p.cX)

    System(X, rr, dXdr, hX, z, hz )
end

# 2nd order accurate finite differences for 1st order derivatives along a
# periodic direction
function Dz!(f_z, f, sys::System)
    NX, Nz = size(f)
    odz2 = 0.5 / sys.hz

    @inbounds for j in 2:Nz-1
        @inbounds for i in 1:NX
            f_z[i,j] = (f[i,j+1] - f[i,j-1]) * odz2
        end
    end

    @inbounds for i in 1:NX
        f_z[i,1]   = (f[i,2] - f[i,end]) * odz2
        f_z[i,end] = (f[i,1] - f[i,end-1]) * odz2
    end

    f_z
end
function Dz(f, sys)
    f_z = similar(f)
    Dz!(f_z, f, sys)
end


# 2nd order accurate finite difference operator for 1st order derivatives along
# a non-periodic direction. forward and backward 2nd order accurate finite
# differences are used on the first and last points of the grid respectively
function DX!(f_X, f, sys::System)
    NX, Nz = size(f)
    odX2 = 0.5 / sys.hX

    @inbounds for j in 1:Nz
        @inbounds for i in 2:NX-1
            f_X[i,j] = (f[i+1,j] - f[i-1,j]) * odX2
        end
    end

    @inbounds for j in 1:Nz
        f_X[1,j]   = (-3.0* f[1,j]   + 4.0*f[2,j]     - f[3,j]) * odX2
        f_X[end,j] = ( 3.0* f[end,j] - 4.0*f[end-1,j] + f[end-2,j]) * odX2
    end

    f_X
end
function DX(f, sys)
    f_X = similar(f)
    DX!(f_X, f, sys)
end

# this defines the functions that determine the BCs of ϕ and ψv at r=rmin
function setup_BC(p::Param)
    amp = p.noise_amplitude
    global ϕ0_of_uz = (u,z) -> randn(Float64)*(amp) #*sin(z)
    global ψv0_of_uz = (u,z) -> randn(Float64)*(amp) #*sin(z)
end

# Initial data: the ingoing mode at initial time i.e. ψ(u=0)
function init_ψ(sys::System, p::Param)
    NX = length(sys.X)
    Nz = length(sys.z)

    amp = p.noise_amplitude
    randn(NX, Nz) * amp
end

# integrate in the null hypersurface
function get_ϕψv!(ϕ, ψv, ψ, u, sys)
    NX, Nz = size(ψ)

    S_ϕ  = copy(ψ)

    Threads.@threads for j in 1:Nz
        itp  = interpolate( S_ϕ[:,j], BSpline(Cubic(Flat(OnGrid()))) )
        sitp = scale(itp, sys.X[1]:sys.hX:sys.X[end])
        rhs_ϕ!(f, p, x) = sitp(x)

        Xspan = (sys.X[1], sys.X[end])

        # this defines the outgoing mode i.e. boundary condition at each
        # timestep
        ϕ0 = ϕ0_of_uz(u,sys.z[j])
        # define the PDE problem
        prob_ϕ = ODEProblem(rhs_ϕ!, ϕ0, Xspan)
        # solve the PDE problem
        sol_ϕ  = solve(prob_ϕ, SSPRK22(), dt=sys.hX, adaptive=false)
        # pass the above solution
        ϕ[:, j] .= sol_ϕ.(sys.X)
    end

    S_ψv  = ϕ .+ ψ

    Threads.@threads for j in 1:Nz
        itp  = interpolate( S_ψv[:,j], BSpline(Cubic(Flat(OnGrid()))) )
        sitp = scale(itp, sys.X[1]:sys.hX:sys.X[end])
        rhs_ψv!(f, p, x) = sitp(x)

        Xspan = (sys.X[1], sys.X[end])

        # this defines the outgoing mode i.e. boundary condition at each
        # timestep
        ψv0 = ψv0_of_uz(u,sys.z[j])
        # define the PDE problem
        prob_ψv = ODEProblem(rhs_ψv!, ψv0, Xspan)
        # solve the PDE problem
        sol_ψv  = solve(prob_ψv, SSPRK22(), dt=sys.hX, adaptive=false)
        # pass the above solution
        ψv[:, j] .= sol_ψv.(sys.X)
    end

    ϕ, ψv
end
function get_ϕψv(ψ, u, sys)
    ϕ    = similar(ψ)
    ψv   = similar(ψ)
    get_ϕψv!(ϕ, ψv, ψ, u, sys)
end


# the rhs for the time evolution equation
function get_ψ_u!(ψ_u, ψ, ϕ, sys)
    ψ_X = DX(ψ, sys)
    ψ_z = Dz(ψ, sys)
    ψ_r = sys.dXdr .* ψ_X

    ψ_u .= 0.5 * ψ_r .+ ψ_z .+ ϕ
end
function rhs_ψ!(dψ, ψ, sys, u)
    ϕ, ψv = get_ϕψv(ψ, u, sys)
    get_ψ_u!(dψ, ψ, ϕ, sys)
    nothing
end

function write_2D(it::Int, t, data_dir::String, ψ::Array, ψv::Array, ϕ::Array)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "data_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "ψ",  ψ)
        write(file, "ψv", ψv)
        write(file, "ϕ",  ϕ)
        attrs(file)["time"] = t
    end
    nothing
end

# function that performs the time evolution
function run_toy(p::Param)

    # pass the parameters of the system
    sys = System(p)

    # create the folders where data are saved
    data_dir = mkpath(p.out_dir)

    # setup BC functions
    setup_BC(p)

    # initialize ψ
    ψ = init_ψ(sys, p)

    # time span of the simulation
    tspan = (0.0, p.umax)

    # timestep. careful with CFL condition
    dt0   = 0.25 * minimum([sys.hX, sys.hz])

    # define the PDE problem for time integration
    prob  = ODEProblem(rhs_ψ!, ψ, tspan, sys)
    # http://docs.juliadiffeq.org/latest/basics/integrator.html
    integrator = init(prob, RK4(), save_everystep=false, dt=dt0, adaptive=false)

    # write the coordinates. if there's a problem with this for windows users,
    # it can be done differently (with joinpath, or similar)
    h5write(data_dir*"/x.h5", "x", sys.X)
    h5write(data_dir*"/z.h5", "z", sys.z)
    h5write(data_dir*"/r.h5", "r", sys.r)

    it = 0
    t  = 0.0
    ϕ, ψv = get_ϕψv(ψ, t, sys)

    # save initial data
    write_2D(it, t, data_dir, ψ, ψv, ϕ)

    if p.compute_L2_norm
        L2_norm   = sqrt(sys.hz*sys.hX*sum(ψ.*ψ  + ψv.*ψv + ϕ.*ϕ))
        L2_norm_t = [[t, L2_norm]]
    end

    if p.compute_Lop_norm
        ϕz = Dz(ϕ, sys)
        Lopside_norm   = sqrt(sys.hz*sys.hX*sum(ψ.*ψ  + ψv.*ψv + ϕ.*ϕ + ϕz.*ϕz ))
        Lopside_norm_t = [[t, Lopside_norm]]
    end

    println("-------------------------------------------------------------")
    println("Iteration      Time |            ψ ")
    println("                    |    minimum      maximum")
    println("-------------------------------------------------------------")
    @printf "%9d %9.3f |  %9.4g    %9.4g\n" it t minimum(ψ) maximum(ψ)

    # start time evolution
    for (f,t) in tuples(integrator)
        it += 1

        ψ = f
        get_ϕψv!(ϕ, ψv, ψ, t, sys)

        @printf "%9d %9.3f |  %9.4g    %9.4g\n" it t minimum(ψ) maximum(ψ)

        if p.compute_L2_norm
            L2_norm   = sqrt(sys.hz*sys.hX*sum(ψ.*ψ  + ψv.*ψv + ϕ.*ϕ))
            push!(L2_norm_t, [t, L2_norm])
        end

        if p.compute_Lop_norm
            Dz!(ϕz, ϕ, sys)
            Lopside_norm   = sqrt(sys.hz*sys.hX*sum(ψ.*ψ  + ψv.*ψv + ϕ.*ϕ + ϕz.*ϕz ))
            push!(Lopside_norm_t, [t, Lopside_norm])
        end

        if it % p.out_every == 0
            write_2D(it, t, data_dir, ψ, ψv, ϕ)
        end

    end

    # store norms

    if p.compute_L2_norm
        outfile = joinpath(data_dir, "L2_norm.dat")
        open(outfile, "w") do io
            println(io, "# t      |      norm_L2")
            writedlm(io, L2_norm_t)
        end
    end

    if p.compute_Lop_norm
        outfile = joinpath(data_dir, "Lopside_norm.dat")
        open(outfile, "w") do io
            println(io, "# t      |      norm_Lop")
            writedlm(io, Lopside_norm_t)
        end
    end

    println("-------------------------------------------------------------")
    println("Done.")

end

end # module
