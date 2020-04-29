
using BondiToy
using Parameters

struct SH_smooth_B1 <: NestedIBVP end

# TODO: maybe promote 0.5 and 0.1 to parameters?
function outgoing(u, z)
    exp( -((u - 0.5) / 0.1)^2 ) * sin(z)
end


BondiToy.ϕ0_of_uz(u::T, z::T, ibvp::NestedIBVP) where {T<:Real} =
    3 * outgoing(u, z)

BondiToy.ψv0_of_uz(u::T, z::T, ibvp::NestedIBVP) where {T<:Real} =
    outgoing(u, z)

# smooth initial data: gaussian in X and sin in z
# TODO: maybe promote 0.5 and 0.1 to parameters?
BondiToy.ψ0_of_Xz(X::T, z::T, ibvp::NestedIBVP) where {T<:Real} =
    exp( -((X - 0.5)/0.1)^2 ) * sin(z)


toy_model = "SH_smooth_B1"
root_dir="./run00"

D = 1
noise_amplitude_drop = 0.25

NX = (17-1)*2^D + 1 #17 coarse
Nz = 16*2^D #16 coarse

# parameters to be passed in the toy model
p = Param(
    NX = NX,
    Nz = Nz,
    cX = 10.0,
    rmin = 2.0,
    umax = 1.0, # total time of the simulation in code units

    # uncomment the appropriate SCALING AMPLITUDE for noisy given data
    # for L2 norm
    # for the LOPSIDED norm
    #noise_amplitude_drop = 0.125,

    out_dir = joinpath(root_dir, toy_model*"_"*string(noise_amplitude_drop),
                       "data_$(NX)_$(Nz)"),

    out_every = 4*16*2^D,
)

ibvp = SH_smooth_B1()

run_toy(p, ibvp)
