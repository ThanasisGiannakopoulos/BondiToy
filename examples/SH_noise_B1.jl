
using BondiToy

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

    noise_amplitude = noise_amplitude_drop^D,

    out_dir = joinpath(root_dir, toy_model*"_"*string(noise_amplitude_drop),
                       "data_$(NX)_$(Nz)"),

    out_every = 4*16*2^D,
)

run_toy(p)
