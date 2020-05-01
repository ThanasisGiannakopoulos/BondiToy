
using BondiToy
using Parameters
using Random

@with_kw struct WH_noise_B1 <: NestedIBVP
    noise_amplitude :: Float64
    # the element of the angular principal matrix Az that controls the
    # hyperbolicity of the toy model.
    az_21 :: Int # 0 for SH; 1 for WH
    # the elements of the matrix B that encodes the sourc terms S,
    # with S = B*state_vector
    b11 :: Int
    b13 :: Int
    b21 :: Int
    b22 :: Int
    b23 :: Int
    b31 :: Int
    b32 :: Int
    b33 :: Int
end

BondiToy.ϕ0_of_uz(u::T, z::T, ibvp::WH_noise_B1) where {T<:Real} =
    ibvp.noise_amplitude * randn(T)

BondiToy.ψv0_of_uz(u::T, z::T, ibvp::WH_noise_B1) where {T<:Real} =
    ibvp.noise_amplitude * randn(T)

BondiToy.ψ0_of_Xz(X::T, z::T, ibvp::WH_noise_B1) where {T<:Real} =
    ibvp.noise_amplitude * randn(T)


toy_model = "WH_noise_B1_norms"
root_dir="./run00"

D = 4
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

ibvp = WH_noise_B1(
    noise_amplitude = noise_amplitude_drop^D,
    az_21 = 1, # 0 for SH; 1 for WH
    b11 = 0,
    b13 = 1,
    b21 = 1,
    b22 = 0,
    b23 = 1,
    b31 = 1,
    b32 = 0,
    b33 = 0,
)

run_toy(p, ibvp)
