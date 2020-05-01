using HDF5
using DelimitedFiles

# CHANGE THIS to the directory with the data of the model for which
# you want to compute the L2 self norms i.e. the ones used in the self
# convergence ratio
#data_dir = "/home/thanasis/repos/BondiToy/examples/run00/WH_smooth_B0"
data_dir = "/home/thanasis/repos/BondiToy/examples/run00/SH_smooth_B1"

mkpath(data_dir*"/norms_self_2")
data_dir2 = data_dir*"/norms_self_2"
Nmax = 5 # CHANGE SUITABLY  the maximum times we double resolution

Nx = 17 # the overal course graining
Nz = 16

for n in 0:1:Nmax-2 # we need 3 different resolutions to buld the L2
                    # that is used in the self convergence ratio


    # load the x grid
    xc = h5read(data_dir*"/data_$((Nx-1)*2^n + 1)_$(Nz*2^n)/x.h5","x")
    # define dx of each coarse resolution
    dxc = xc[2] - xc[1] 
    xm = h5read(data_dir*"/data_$((Nx-1)*2^(n+1) + 1)_$(Nz*2^(n+1))/x.h5","x")
    xh = h5read(data_dir*"/data_$((Nx-1)*2^(n+2) + 1)_$(Nz*2^(n+2))/x.h5","x")

    # load the z grid
    zc = h5read(data_dir*"/data_$((Nx-1)*2^n + 1)_$(Nz*2^n)/z.h5","z")
    # define dz of each coarse resolution
    dzc = zc[2] - zc[1]
    zm = h5read(data_dir*"/data_$((Nx-1)*2^(n+1) + 1)_$(Nz*2^(n+1))/z.h5","z")
    zh = h5read(data_dir*"/data_$((Nx-1)*2^(n+2) + 1)_$(Nz*2^(n+2))/z.h5","z")

    # load all the timesteps that data is written for
    tc = readdlm(data_dir*"/data_$((Nx-1)*2^n + 1)_$(Nz*2^n)/L2_norm.dat", comments=true)[:,1]
    tm = readdlm(data_dir*"/data_$((Nx-1)*2^(n+1) + 1)_$(Nz*2^(n+1))/L2_norm.dat",
                 comments=true)[:,1]
    th = readdlm(data_dir*"/data_$((Nx-1)*2^(n+2) + 1)_$(Nz*2^(n+2))/L2_norm.dat",
                 comments=true)[:,1]

    # create lists to save the common timesteps of the 3 resolutions
    tm_c = zeros(length(tc))
    th_c = zeros(length(tc))

    # create the lists for the L2 norms of the differences course-med and med-high
    # empty list, we append after
    L2_cmt = zeros(0) 
    L2_mht = zeros(0)

    for i in 0:1:length(tc)-1

        # indices to take care of common timesteps
        i0 = i
        i1 = 2*i0
        i2 = 2*i1
        it_str_0  = lpad(i0, 4, "0")
        it_str_1  = lpad(i1, 4, "0")
        it_str_2  = lpad(i2, 4, "0")

        tm_c[i0+1] = tm[i1+1] - tc[i0+1]
        th_c[i0+1] = th[i2+1] - tc[i0+1]

        println("t = $(tc[i0+1])")
        
        ψc =  h5read(data_dir*"/data_$((Nx-1)*2^n + 1)_$(Nz*2^n)/data_$(it_str_0).h5","ψ")
        ψvc =  h5read(data_dir*"/data_$((Nx-1)*2^n + 1)_$(Nz*2^n)/data_$(it_str_0).h5","ψv")
        ϕc =  h5read(data_dir*"/data_$((Nx-1)*2^n + 1)_$(Nz*2^n)/data_$(it_str_0).h5","ϕ")
        
        ψm =  h5read(data_dir*"/data_$((Nx-1)*2^(n+1) + 1)_$(Nz*2^(n+1))/data_$(it_str_1).h5","ψ")
        ψvm =  h5read(data_dir*"/data_$((Nx-1)*2^(n+1) + 1)_$(Nz*2^(n+1))/data_$(it_str_1).h5","ψv")
        ϕm =  h5read(data_dir*"/data_$((Nx-1)*2^(n+1) + 1)_$(Nz*2^(n+1))/data_$(it_str_1).h5","ϕ")

        ψh =  h5read(data_dir*"/data_$((Nx-1)*2^(n+2) + 1)_$(Nz*2^(n+2))/data_$(it_str_2).h5","ψ")
        ψvh =  h5read(data_dir*"/data_$((Nx-1)*2^(n+2) + 1)_$(Nz*2^(n+2))/data_$(it_str_2).h5","ψv")
        ϕh =  h5read(data_dir*"/data_$((Nx-1)*2^(n+2) + 1)_$(Nz*2^(n+2))/data_$(it_str_2).h5","ϕ")
        
        ψcm = zeros(0)
        ψmh = zeros(0)
        ψvcm = zeros(0)
        ψvmh = zeros(0)
        ϕcm = zeros(0)
        ϕmh = zeros(0)
        
        for j in 1:1:length(xc)-1

            # indices to take care of common x-grid points
            l0 = j
            l1 = 2*(l0-1) + 1
            l2 = 2*(l1-1) + 1

            for p in 1:1:length(zc)

                # indiced to take care of common z-grid points
                p0 = p
                p1 = 2*p0 - 1
                p2 = 2*p1 - 1
            
                append!(ψcm, ψc[l0, p0] - ψm[l1, p1])
                append!(ψmh, ψm[l1, p1] - ψh[l2, p2])
                append!(ψvcm, ψvc[l0, p0] - ψvm[l1, p1])
                append!(ψvmh, ψvm[l1, p1] - ψvh[l2, p2])
                append!(ϕcm, ϕc[l0, p0] - ϕm[l1, p1])
                append!(ϕmh, ϕm[l1, p1] - ϕh[l2, p2])

            end
            
        end

        for p in 1:1:length(zc)
            p0 = p
            p1 = 2*p0 - 1
            p2 = 2*p1 - 1
            
            append!(ψcm, ψc[end, p0] - ψm[end, p1])
            append!(ψmh, ψm[end, p1] - ψh[end, p2])
            append!(ψvcm, ψvc[end, p0] - ψvm[end, p1])
            append!(ψvmh, ψvm[end, p1] - ψvh[end, p2])
            append!(ϕcm, ϕc[end, p0] - ϕm[end, p1])
            append!(ϕmh, ϕm[end, p1] - ϕh[end, p2])

        end
        
        append!( L2_cmt, sqrt(sum(dzc*dxc*( ψcm.*ψcm + ψvcm.*ψvcm + ϕcm.*ϕcm))) )
        append!( L2_mht, sqrt(sum(dzc*dxc*( ψmh.*ψmh + ψvmh.*ψvmh + ϕmh.*ϕmh))) )

    end

    h5write(data_dir2*"/L2_$(n)_$(n+1)_c$(n).h5","L2",L2_cmt)
    h5write(data_dir2*"/L2_$(n+1)_$(n+2)_c$(n).h5","L2",L2_mht)

    h5write(data_dir2*"/t_$(n).h5","t",tc)
    h5write(data_dir2*"/t_$(n)_$(n+1).h5","t",tm_c)
    h5write(data_dir2*"/t_$(n)_$(n+2).h5","t",th_c)

end
