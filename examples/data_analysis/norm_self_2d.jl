
using HDF5
using DelimitedFiles

function list_h5_files(foldername::String; prefix::String="data_")
    path     = abspath(foldername)
    allfiles = readdir(path)

    Ns = length(prefix)

    its_names = Tuple[]
    # append only the files whose names start with the given prefix
    for file in allfiles
        try
            if (file[1:Ns] == prefix && file[end-2:end] == ".h5")
                # extract iteration
                it_str = file[Ns+1:end-3]
                fullname = joinpath(path, file)
                # add to list of tuples with iteration and name
                push!(its_names, (parse(Int64, it_str), fullname))
            end
        catch ex
            if isa(ex, BoundsError)
                # probably triggered by string comparison; do nothing
            else
                throw(ex)
            end
        end
    end

    # sort according to iteration
    sort!(its_names)
    # and extract the list of filenames and iterations
    filenames = [name for (it, name) in its_names]
    its       = [it for (it, name) in its_names]

    (its, filenames)
end


function L2_cmh_t(dir_c, dir_m, dir_h)
    # x grid
    xc = h5read(dir_c * "/x.h5", "x")
    xm = h5read(dir_m * "/x.h5", "x")
    xh = h5read(dir_h * "/x.h5", "x")
    # z grid
    zc = h5read(dir_c * "/z.h5", "z")
    zm = h5read(dir_m * "/z.h5", "z")
    zh = h5read(dir_h * "/z.h5", "z")

    # make sure that we can inject points from the medium and high resolution
    # grids in the coarse grid without interpolation
    @assert xc ≈ xm[1:2:end] ≈ xh[1:4:end]
    @assert zc ≈ zm[1:2:end] ≈ zh[1:4:end]

    # dx of coarse resolution
    dxc = xc[2] - xc[1] 
    # dz of coarse resolution
    dzc = zc[2] - zc[1]

    # list all available iterations (and corresponding files)
    (its_c, all_filenames_c) = list_h5_files(dir_c, prefix="data_")
    (its_m, all_filenames_m) = list_h5_files(dir_m, prefix="data_")
    (its_h, all_filenames_h) = list_h5_files(dir_h, prefix="data_")

    # we want to compare the common timesteps. we assume here that there is a
    # factor of 2 between the lowest resolution and the corresponding higher
    # resolution ones
    filenames_c = all_filenames_c[:]
    filenames_m = all_filenames_m[1:2:end]
    filenames_h = all_filenames_h[1:4:end]

    Nf = length(filenames_c)
    @assert length(filenames_m) == length(filenames_h) == Nf

    tt     = zeros(Nf)
    L2_cmt = zeros(Nf)
    L2_mht = zeros(Nf)

    # initiate the grid function for the outgoing norm
    l2_out_cmt_x = zeros(length(xc))
    l2_out_mht_x = zeros(length(xc))
    
    for it in 1:Nf
        file_c = filenames_c[it]
        file_m = filenames_m[it]
        file_h = filenames_h[it]

        ψc  = h5read(file_c, "ψ")
        ψvc = h5read(file_c, "ψv")
        ϕc  = h5read(file_c, "ϕ")
        tc  = h5readattr(file_c, "./")["time"]

        ψm  = h5read(file_m, "ψ")
        ψvm = h5read(file_m, "ψv")
        ϕm  = h5read(file_m, "ϕ")
        tm  = h5readattr(file_m, "./")["time"]

        ψh  = h5read(file_h, "ψ")
        ψvh = h5read(file_h, "ψv")
        ϕh  = h5read(file_h, "ϕ")
        th  = h5readattr(file_h, "./")["time"]

        # make sure we're comparing the same timestep
        @assert tc ≈ tm ≈ th

        # compute the differences between resolutions (projected onto the coarsest grid)
        ψcm  = ψc  -  ψm[1:2:end, 1:2:end]
        ψvcm = ψvc - ψvm[1:2:end, 1:2:end]
        ϕcm  = ϕc  -  ϕm[1:2:end, 1:2:end]
        ψmh  =  ψm[1:2:end, 1:2:end]  - ψh[1:4:end, 1:4:end]
        ψvmh = ψvm[1:2:end, 1:2:end] - ψvh[1:4:end, 1:4:end]
        ϕmh  =  ϕm[1:2:end, 1:2:end]  - ϕh[1:4:end, 1:4:end]

        tt[it]     = tc

        # define the timestep; needed for the sum in u
        dt0 = 0.25*minimum([dxc, dzc])

        # get the outgoing part of the norm
        l2_out_cmt_x += dt0*dzc*sum( ψvcm.*ψvcm + ϕcm.*ϕcm, dims=2)
        l2_out_mht_x += dt0*dzc*sum( ψvmh.*ψvmh + ϕmh.*ϕmh, dims=2)

        # sum the outgoing and ingoing norms to get the complete one
        L2_cmt[it] = sqrt( sum(dzc*dxc*( ψcm.*ψcm ))) + maximum(sqrt.(l2_out_cmt_x))
        L2_mht[it] = sqrt( sum(dzc*dxc*( ψmh.*ψmh ))) + maximum(sqrt.(l2_out_mht_x))
    end

    tt, L2_cmt, L2_mht
end


# maximum times we double resolution
Nmax = 5

# base lowest resolution
Nx = 17
Nz = 16

root_dir  = "/home/thanasis/repos/BondiToy/examples/run00/"
#"/home/mzilhao/dev/julia/BondiToy/examples/run00/"
toy_model = "SH_smooth_B1/"
    #"WH_smooth_B0/"

# we need 3 different resolutions to build the L2 norm that is used in the self
# convergence ratio
for n in 0:1:Nmax-2
    dir_c = joinpath(root_dir, toy_model, "data_$((Nx-1)*2^n + 1)_$(Nz*2^n)")
    dir_m = joinpath(root_dir, toy_model, "data_$((Nx-1)*2^(n+1) + 1)_$(Nz*2^(n+1))")
    dir_h = joinpath(root_dir, toy_model, "data_$((Nx-1)*2^(n+2) + 1)_$(Nz*2^(n+2))")

    tt, L2_cmt, L2_mht = L2_cmh_t(dir_c, dir_m, dir_h)

    data_dir2 = joinpath(root_dir, toy_model, "norms_self_2")
    mkpath(data_dir2)

    outfile  = joinpath(data_dir2, "L2_$(n)_$(n+1)_c$(n).dat")
    open(outfile, "w") do io
        println(io, "# t      |      L2")
        writedlm(io, [tt L2_cmt])
    end

    outfile  = joinpath(data_dir2, "L2_$(n+1)_$(n+2)_c$(n).dat")
    open(outfile, "w") do io
        println(io, "# t      |      L2")
        writedlm(io, [tt L2_mht])
    end
end
