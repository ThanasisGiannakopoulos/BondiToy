module BondiToy

export Param, NestedIBVP, CoupledIBVP, run_toy

using Parameters
using DifferentialEquations
using Interpolations
using HDF5
using Printf
using RecursiveArrayTools

@with_kw struct Param
    # discretization parameters
    NX                :: Int
    Nz                :: Int
    cX                :: Float64 = 10.0
    rmin              :: Float64 = 2.0
    umax              :: Float64
    # directory to save data
    out_dir           :: String

    # how often to save data
    out_every               :: Int
    compute_L2_norm_every   :: Int = 1 # set to negative values to turn off
    compute_Lop_norm_every  :: Int = 1 # set to negative values to turn off
end


"""
Extend this type for the initial and boundary condition parameters. The IBVP can
be nested (each intrinsic equation can be solved one after the other) or coupled
(the intrinsic equations are coupled).
"""
abstract type IBVP end

"""
The NestedIBVP structure encodes all the models with nested intrinsic
equations, both SH and WH, and for sources B0, B1 and B2.  It needs
the paremeters az_21, b11, b13, b21, b22, b33, which are the elements
of the angular principal matrix Az and the lower order matrix B,
respectively.
"""
abstract type NestedIBVP  <: IBVP end
abstract type CoupledIBVP <: IBVP end

"""
Function that defines the BCs of ϕ at r=rmin

Needs signature ϕ0_of_uz(u::T, z::T, ibvp::IBVP) where {T<:Real}
"""
function ϕ0_of_uz end


"""
Function that defines the BCs of ψv at r=rmin

Needs signature ψv0_of_uz(u::T, z::T, ibvp::IBVP) where {T<:Real}
"""
function ψv0_of_uz end


"""
Function that defines the ingoing mode of ψ at u=0, i.e. ψ(u=0,X,z)

Needs signature ψ0_of_Xz(X::T, z::T, ibvp::IBVP) where {T<:Real}
"""
function ψ0_of_Xz end


rtoX(r, cX, rmin) = (r .- rmin) ./ sqrt.(cX*cX .+ (r .- rmin).^2)
Xtor(X, cX, rmin) = rmin .+ cX * X ./ sqrt.(1.0 .- X.^2)
dX_dr(X, cX) = (1.0 .- X.^2).^1.5 ./ cX

struct System{T}
    X    :: Vector{T}
    r    :: Vector{T}
    dXdr :: Vector{T}
    hX   :: T
    z    :: Vector{T}
    hz   :: T
end
function System(p::Param)
    hX = 1.0 / (p.NX-1)
    X  = range(0, 1, length=p.NX)

    hz = 2.0*π / (p.Nz) # remove last point; 0=2π (periodic)
    z  = range(0, 2.0*π-hz, length=p.Nz)

    rr    = Xtor(X, p.cX, p.rmin)
    dXdr  = dX_dr(X, p.cX)

    System{typeof(hX)}(X, rr, dXdr, hX, z, hz)
end

# 2nd order accurate finite differences for 1st order derivatives along a
# periodic direction
function Dz!(f_z::Matrix, f::Matrix, sys::System)
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

# now acting on vectors
function Dz!(f_z::Vector, f::Vector, sys::System)
    Nz = length(f)
    odz2 = 0.5 / sys.hz

    @inbounds for j in 2:Nz-1
        f_z[j] = (f[j+1] - f[j-1]) * odz2
    end

    f_z[1]   = (f[2] - f[end]) * odz2
    f_z[end] = (f[1] - f[end-1]) * odz2

    f_z
end

function Dz(f, sys::System)
    f_z = similar(f)
    Dz!(f_z, f, sys)
end

# 2nd order accurate finite difference operator for 1st order derivatives along
# a non-periodic direction. forward and backward 2nd order accurate finite
# differences are used on the first and last points of the grid respectively
function DX!(f_X::Matrix, f::Matrix, sys::System)
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
function DX(f, sys::System)
    f_X = similar(f)
    DX!(f_X, f, sys)
end

function init_ψ(sys::System{T}, ibvp::IBVP) where {T}
    NX = length(sys.X)
    Nz = length(sys.z)

    ψ0 = zeros(T, NX, Nz)

    @inbounds for j in 1:Nz
        @inbounds for i in 1:NX
            ψ0[i,j] = ψ0_of_Xz(sys.X[i], sys.z[j], ibvp)
        end
    end
    ψ0
end

# integrate in the null hypersurface for nested intrinsic eqs
function get_ϕψv!(ϕ, ψv, ψ, u, sys::System, ibvp::NestedIBVP)
    NX, Nz = size(ψ)

    S_ϕ  = ibvp.b13 * ψ

    Threads.@threads for j in 1:Nz
        itp  = interpolate( S_ϕ[:,j], BSpline(Cubic(Flat(OnGrid()))) )
        sitp = scale(itp, sys.X[1]:sys.hX:sys.X[end])
        rhs_ϕ(f, p, x) = sitp(x) + ibvp.b11 * f

        Xspan = (sys.X[1], sys.X[end])

        # this defines the outgoing mode i.e. boundary condition at each
        # timestep
        ϕ0 = ϕ0_of_uz(u, sys.z[j], ibvp)
        # define the PDE problem
        prob_ϕ = ODEProblem(rhs_ϕ, ϕ0, Xspan)
        # solve the PDE problem
        sol_ϕ  = solve(prob_ϕ, SSPRK22(), dt=sys.hX, adaptive=false)
        # pass the above solution
        ϕ[:, j] = sol_ϕ.(sys.X)
    end

    ϕ_z = Dz(ϕ, sys)
    S_ψv  = ibvp.az_21 * ϕ_z .+ ibvp.b21 * ϕ .+ ibvp.b23 * ψ

    Threads.@threads for j in 1:Nz
        itp  = interpolate( S_ψv[:,j], BSpline(Cubic(Flat(OnGrid()))) )
        sitp = scale(itp, sys.X[1]:sys.hX:sys.X[end])
        rhs_ψv(f, p, x) = sitp(x) + ibvp.b22 * f

        Xspan = (sys.X[1], sys.X[end])

        # this defines the outgoing mode i.e. boundary condition at each
        # timestep
        ψv0 = ψv0_of_uz(u, sys.z[j], ibvp)
        # define the PDE problem
        prob_ψv = ODEProblem(rhs_ψv, ψv0, Xspan)
        # solve the PDE problem
        sol_ψv  = solve(prob_ψv, SSPRK22(), dt=sys.hX, adaptive=false)
        # pass the above solution
        ψv[:, j] = sol_ψv.(sys.X)
    end

    ϕ, ψv
end

# integrate in the null hypersurface for non-nested intrinsic eqs
function intrinsic_rhs!(dv, v, sys, x)
    dϕ  = dv[1]
    dψv = dv[2]
    ϕ   = v[1]
    ψv  = v[2]

    # the rhs for the intrinsic coupled PDE of ϕ and ψv
    dϕ .= ψv
    Dz!(dψv, ϕ, sys)  # dψv = ∂_z ϕ
    dv
end
function get_ϕψv!(ϕ, ψv, ψ, u, sys::System{T}, ibvp::CoupledIBVP) where {T}
    NX, Nz = size(ψ)
    ϕ0  = zeros(T, Nz)
    ψv0 = zeros(T, Nz)

    # this defines the outgoing mode i.e. boundary condition at each timestep
    @inbounds for j in 1:Nz
        ϕ0[j]  = ϕ0_of_uz(u, sys.z[j], ibvp)
        ψv0[j] = ψv0_of_uz(u, sys.z[j], ibvp)
    end

    # save the outgoing mode at X=X0 to ϕ and ψv of (x,z)
    ϕ[1,:]  = ϕ0
    ψv[1,:] = ψv0

    # write the ID as a vector for the coupled PDE system
    v0 = VectorOfArray([ϕ0, ψv0])

    # define the x-span to solve the PDE for
    Xspan = (sys.X[1], sys.X[end])

    # write the PDE as an ODE prob
    intrinsic_prob = ODEProblem(intrinsic_rhs!, v0, Xspan, sys)

    # def the x-integrator
    x_integrator = init(intrinsic_prob, SSPRK22(),
                        save_everytimestep=false, dt=sys.hX, adaptive=false)

    # the following assumes that the integrator is stepping *exactly* at the X
    # grid points. this should be true with the above options of dt=sys.hX,
    # adaptive=false, but we must remember to change things if using an adaptive
    # integrator or something similar.
    iter = 1
    for (f,t) in tuples(x_integrator)
        iter += 1

        ϕ[iter,:]  = f[1]
        ψv[iter,:] = f[2]
    end

    ϕ, ψv
end
function get_ϕψv(ψ, u, sys, ibvp)
    ϕ    = similar(ψ)
    ψv   = similar(ψ)
    get_ϕψv!(ϕ, ψv, ψ, u, sys, ibvp)
end


# the rhs for the time evolution equation
function get_ψ_u!(ψ_u, ψ, ϕ, ψv, sys, ibvp)
    ψ_X = DX(ψ, sys)
    ψ_z = Dz(ψ, sys)
    ψ_r = sys.dXdr .* ψ_X

    ψ_u .= 0.5 * ψ_r .+ ψ_z .+ ibvp.b31 * ϕ .+ ibvp.b32 * ψv .+ ibvp.b33 * ψ
end
function rhs_ψ!(dψ, ψ, (sys, ibvp), u)
    ϕ, ψv = get_ϕψv(ψ, u, sys, ibvp)
    get_ψ_u!(dψ, ψ, ϕ, ψv, sys, ibvp)
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
function run_toy(p::Param, ibvp::IBVP)

    # pass the parameters of the system
    sys = System(p)

    # create the folders where data are saved
    data_dir = mkpath(p.out_dir)

    # initialize ψ
    ψ = init_ψ(sys, ibvp)

    # time span of the simulation
    tspan = (0.0, p.umax)

    # timestep. careful with CFL condition
    dt0   = 0.25 * minimum([sys.hX, sys.hz])

    # define the PDE problem for time integration
    prob  = ODEProblem(rhs_ψ!, ψ, tspan, (sys, ibvp))
    # http://docs.juliadiffeq.org/latest/basics/integrator.html
    integrator = init(prob, RK4(), save_everystep=false, dt=dt0, adaptive=false)

    # write the coordinates. if there's a problem with this for windows users,
    # it can be done differently (with joinpath, or similar)
    h5write(data_dir*"/x.h5", "x", sys.X)
    h5write(data_dir*"/z.h5", "z", sys.z)
    h5write(data_dir*"/r.h5", "r", sys.r)

    it = 0
    t  = 0.0
    ϕ, ψv = get_ϕψv(ψ, t, sys, ibvp)

    # save initial data
    write_2D(it, t, data_dir, ψ, ψv, ϕ)

    
    if p.compute_L2_norm_every > 0
        # initiate the grid function to be used for the sum in the time domain
        # here we sum over the z domain for the initial time
        l2_out_x  = dt0*sys.hz*sum(ψv.*ψv + ϕ.*ϕ, dims=2)
        L2_norm = sqrt(sys.hz*sys.hX*sum(ψ.*ψ)) + maximum(sqrt.(l2_out_x))
        outfile = joinpath(data_dir, "L2_norm.dat")
        open(outfile, "w") do io
            println(io, "# t        |      norm_L2")
            println(io, "$t  \t  $L2_norm")
        end
    end

    if p.compute_Lop_norm_every > 0
        ϕz = Dz(ϕ, sys)
        # initiate the grid function to be used for the sum in the time domain
        # here we sum over the z domain for the initial time
        lop_out_x = dt0*sys.hz*sum(ψv.*ψv + ϕ.*ϕ + ϕz.*ϕz, dims=2)
        Lopside_norm = sqrt(sys.hz*sys.hX*sum(ψ.*ψ)) + maximum(sqrt.(lop_out_x))
        outfile = joinpath(data_dir, "Lopside_norm.dat")
        open(outfile, "w") do io
            println(io, "# t        |      norm_Lop")
            println(io, "$t  \t  $Lopside_norm")
        end
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
        get_ϕψv!(ϕ, ψv, ψ, t, sys, ibvp)

        @printf "%9d %9.3f |  %9.4g    %9.4g\n" it t minimum(ψ) maximum(ψ)

        if p.compute_L2_norm_every > 0 && it % p.compute_L2_norm_every == 0
            # add on the l2_out_x grid function the z sums for this time instant
            l2_out_x + = dt0*sys.hz*sum(ψv.*ψv + ϕ.*ϕ, dims=2)
            L2_norm = sqrt(sys.hz*sys.hX*sum(ψ.*ψ)) + maximum(sqrt.(l2_out_x))
            outfile = joinpath(data_dir, "L2_norm.dat")
            open(outfile, "a") do io
                println(io, "$t  \t  $L2_norm")
            end
        end

        if p.compute_Lop_norm_every > 0 && it % p.compute_Lop_norm_every == 0
            Dz!(ϕz, ϕ, sys)
            # add on the lop_out_x grid function the z sums for this time instant
            lop_out_x += dt0*sys.hz*sum(ψv.*ψv + ϕ.*ϕ + ϕz.*ϕz, dims=2)
            Lopside_norm = sqrt(sys.hz*sys.hX*sum(ψ.*ψ)) + maximum(sqrt.(lop_out_x))
            outfile = joinpath(data_dir, "Lopside_norm.dat")
            open(outfile, "a") do io
                println(io, "$t  \t  $Lopside_norm")
            end
        end

        if it % p.out_every == 0
            write_2D(it, t, data_dir, ψ, ψv, ϕ)
        end
    end

    println("-------------------------------------------------------------")
    println("Done.")
end

end # module
