
using BondiToy
using Parameters
using Random

@with_kw struct WH_B3 <: CoupledIBVP
    noise_amplitude :: Float64
    b31 :: Int
    b32 :: Int
    b33 :: Int
end

BondiToy.ϕ0_of_uz(u::T, z::T, ibvp::WH_B3) where {T<:Real} =
    ibvp.noise_amplitude * randn(T)

BondiToy.ψv0_of_uz(u::T, z::T, ibvp::WH_B3) where {T<:Real} =
    ibvp.noise_amplitude * randn(T)

BondiToy.ψ0_of_Xz(X::T, z::T, ibvp::WH_B3) where {T<:Real} =
    ibvp.noise_amplitude * randn(T)


toy_model = "WH_noise_B3_norms"
root_dir="./run00"

D = 4
noise_amplitude_drop = 0.125 #0.25

NX = (17-1)*2^D + 1 #17 coarse
Nz = 16*2^D #16 coarse

# parameters to be passed in the toy model
p = Param(
    NX = NX,
    Nz = Nz,
    cX = 10.0,
    rmin = 2.0,
    umax = 14.0, # total time of the simulation in code units

    # uncomment the appropriate SCALING AMPLITUDE for noisy given data
    # for L2 norm
    # for the LOPSIDED norm
    #noise_amplitude_drop = 0.125,

    out_dir = joinpath(root_dir, toy_model*"_"*string(noise_amplitude_drop),
                       "data_$(NX)_$(Nz)"),

    out_every = 4*16*2^D,
)

ibvp = WH_B3(
    noise_amplitude = noise_amplitude_drop^D,
    b31 = 0,
    b32 = 0,
    b33 = 0,
)

run_toy(p, ibvp)
