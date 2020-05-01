
using BondiToy
using Parameters

@with_kw struct WH_smooth_B0 <: NestedIBVP
    BD_gauss_mu    :: Float64
    BD_gauss_sigma :: Float64
    ID_gauss_mu    :: Float64
    ID_gauss_sigma :: Float64
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

# smooth boundary data: gaussian in u and sin in z
function outgoing(u, z, ibvp::WH_smooth_B0)
    exp( -((u - ibvp.BD_gauss_mu) / ibvp.BD_gauss_sigma)^2 ) * sin(z)
end


BondiToy.ϕ0_of_uz(u::T, z::T, ibvp::WH_smooth_B0) where {T<:Real} =
    3 * outgoing(u, z, ibvp)

BondiToy.ψv0_of_uz(u::T, z::T, ibvp::WH_smooth_B0) where {T<:Real} =
    outgoing(u, z, ibvp)

# smooth initial data: gaussian in X and sin in z
BondiToy.ψ0_of_Xz(X::T, z::T, ibvp::WH_smooth_B0) where {T<:Real} =
    exp( -((X - ibvp.ID_gauss_mu) / ibvp.ID_gauss_sigma)^2 ) * sin(z)


toy_model = "WH_smooth_B0"
root_dir="./run00"

D = 0

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

    out_dir = joinpath(root_dir, toy_model,
                       "data_$(NX)_$(Nz)"),

    out_every = 1, #4*16*2^D,
)

ibvp = WH_smooth_B0(
    BD_gauss_mu    = 0.5,
    BD_gauss_sigma = 0.1,
    ID_gauss_mu    = 0.5,
    ID_gauss_sigma = 0.1,
    az_21 = 1, # 0 for SH; 1 for WH
    b11 = 0,
    b13 = 0,
    b21 = 0,
    b22 = 0,
    b23 = 0,
    b31 = 0,
    b32 = 0,
    b33 = 0,
)

run_toy(p, ibvp)
