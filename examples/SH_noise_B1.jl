
using BondiToy
using Parameters
using Random

@with_kw struct SH_noise <: IBVP
    noise_amplitude :: Float64
end

BondiToy.ϕ0_of_uz(u::T, z::T, ibvp::IBVP) where {T<:Real} =
    ibvp.noise_amplitude * randn(T)

BondiToy.ψv0_of_uz(u::T, z::T, ibvp::IBVP) where {T<:Real} =
    ibvp.noise_amplitude * randn(T)

BondiToy.ψ0_of_Xz(X::T, z::T, ibvp::IBVP) where {T<:Real} =
    ibvp.noise_amplitude * randn(T)


toy_model = "SH_noise_B1_norms"
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

ibvp = SH_noise(
    noise_amplitude = noise_amplitude_drop^D,
)

run_toy(p, ibvp)
