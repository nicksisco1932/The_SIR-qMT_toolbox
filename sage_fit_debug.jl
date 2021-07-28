

using Pkg	
try
    println("If is first time you ran the code. It will take a minute to precompile.")
    @eval using NIfTI; 
    @eval using LsqFit;
    @eval using Printf
    @eval using ArgParse;
    @eval using Statistics;
    @eval using Optim;
    Pkg.precompile()
catch e
    # not found; install and try loading again
    Pkg.add("NIfTI")
    Pkg.add("LsqFit")
    Pkg.add("Printf")
    Pkg.add("ArgParse")
    @eval using NIfTI; 
    @eval using LsqFit;
    @eval using Printf
    @eval using ArgParse;
    @eval using Statistics;
    @eval using Optim;
end

using NIfTI; 
using LsqFit;
using Printf
using ArgParse;
using Statistics;
using Optim;

include("./utils.jl")


b_fname = "/mnt/c/Users/nicks/Documents/MRI_data/SAGE/nifti/bPT1319001_preb_mask.nii.gz"
base = "/mnt/c/Users/nicks/Documents/MRI_data/SAGE/nifti/"
paths = [ @sprintf("%sPT1319001_TE1_img_w_Skull.nii.gz",base),
                   @sprintf("%sPT1319001_TE2_img_w_Skull.nii.gz",base),
                               @sprintf("%sPT1319001_TE3_img_w_Skull.nii.gz",base), 
                   @sprintf("%sPT1319001_TE4_img_w_Skull.nii.gz",base),
                   @sprintf("%sPT1319001_TE5_img_w_Skull.nii.gz",base)
                   ]

DATA,MASK,nx,ny,nz,ne,nt=load_sage_data(paths,b_fname);
    tot_voxels = nx*ny*nz
    vec_mask = reshape(MASK[:,:,:,1,1],tot_voxels,1);
    vec_data = reshape(DATA,tot_voxels,ne,nt);

    # MASK = Array{Bool}(niread(b_fname))
    # MASK = MASK.raw
    # vec_mask = reshape(MASK,tot_voxels,1)
    

    


    te1=0.00782;
    te2=0.028769;
    te3=0.060674;
    te4=0.081622;
    te5=0.102571;
    echos=Vector{Float64}([te1,te2,te3,te4,te5]);


    TR = 1800/1000;

    X0=[1000,100.0,50,100];
    IND=findall(x->x.>0,vec_mask[:,1]);

    newDATA=zeros(ne,nx,ny,nz,nt);

    function linearize(nx,ny,nz,ne,nt,data,mask)
        
        # newDATA=zeros(ne,nx,ny,nz,nt);

        # for ii in 1:nx
        #     for jj in 1:ny
        #         for kk in 1:nz
        #             for tt in 1:nt
        #                 if mask[ii,jj,kk]
        #                     newDATA[:,ii,jj,kk,tt] = data[ii,jj,kk,:,tt]
        #                 end
        #             end
        #         end
        #     end
        # end


        # STE1=zeros(nx,ny,nz);
        # STE2=similar(STE1);
        # STE5=similar(STE1);
        # begin
        #     for ii in 1:nx
        #         for jj in 1:ny
        #             for kk in 1:nz
        #                 STE1[ii,jj,kk] = mean(data[ii,jj,kk,1,1:15]);
        #                 STE2[ii,jj,kk] = mean(data[ii,jj,kk,2,1:15]);
        #                 STE5[ii,jj,kk] = mean(data[ii,jj,kk,5,1:15]);
        #             end
        #         end
        #     end
        # end

        # STE1_pre = repeat(STE1,1,1,1,150);
        # STE2_pre = repeat(STE1,1,1,1,150);
        # STE5_pre = repeat(STE1,1,1,1,150);

        ste1 = data[:,:,:,1,:];
        ste2 = data[:,:,:,2,:];
        ste5 = data[:,:,:,5,:];


        
        R₂star_log = similar(ste1);
        R₂_log = similar(ste1);

        IND=findall(x->x==true,mask);
        begin
            for ii in 1:nx
                for jj in 1:ny
                    for kk in 1:nz
                        if mask[ii,jj,kk]
                            for bb = 1:nt
                                a = ste1[ii,jj,kk,bb];
                                # b = STE1_pre[ii,jj,kk,bb];
                                c = ste2[ii,jj,kk,bb];
                                # d = STE2_pre[ii,jj,kk,bb];
                                e = ste5[ii,jj,kk,bb];
                                # f = STE5_pre[ii,jj,kk,bb];
                                ste0 = a*(a/c)^(te1/(te2-te1));
                                R₂star_log[ii,jj,kk,bb] = log( (a) / (c)) /(te2-te1);
                                R₂_log[ii,jj,kk,bb] = log( (ste0 / e)) /te5
                            end
                        end
                    end
                end
            end
        end


        vec_R2fit = zeros(tot_voxels,nt);
        vec_R2sfit = zeros(tot_voxels,nt);
        vec_R2_log = reshape(R₂_log,tot_voxels,nt);
        vec_R2s_log = reshape(R₂star_log,tot_voxels,nt);
        # vec_R2fit = zeros(nx*ny*nz*nt,4);
        # vec_R2sfit = zeros(nx*ny*nz*nt,4);
        # vec_R2_log = reshape(R₂_log,nx*ny*nz*nt);
        # vec_R2s_log = reshape(R₂star_log,nx*ny*nz*nt);
        tmpOUT=Array{Float64}(zeros(nx*ny*nz,4,nt));

        return vec_R2_log,vec_R2s_log,vec_R2sfit,vec_R2fit,tmpOUT
    end

    vec_R2_log,vec_R2s_log,vec_R2sfit,vec_R2fit,tmpOUT = linearize(nx,ny,nz,ne,nt,DATA,MASK)

    function sage_ns_temp(echotimes,p)
        x=echotimes;
        TE=x[end];
        R₂star = p[2]
        R₂ = p[3]
        S₀I = p[1]
        S₀II = p[4]

        M = similar(echotimes)
        for k in 1:length(echotimes)
            if x[k]<=TE/2
                M[k] = S₀I * exp(-x[k] * R₂star)
            elseif x[k]>=TE/2
                M[k] = S₀II * exp(-TE * (R₂star - R₂) * exp(-x[k] * (2*R₂-R₂star)))
            end
        end
        return M
    end

    newDATA=zeros(nx,ny,nz,nt,ne);

    for ii in 1:nx
        for jj in 1:ny
            for kk in 1:nz
                for tt in 1:nt
                    if MASK[ii,jj,kk]
                        newDATA[ii,jj,kk,tt,:] = DATA[ii,jj,kk,:,tt]
                    end
                end
            end
        end
    end
    
    vec_data = reshape(newDATA,nx*ny*nz*nt,ne);
    newMASK = repeat(MASK,1,1,1,150);
    vec_mask = reshape(newMASK,nx*ny*nz*nt);
    tot_voxels=nx*ny*nz*nt;
    tmpOUT = Array{Float64}(zeros(tot_voxels,4));
    x0 = Vector{Float64}(zeros(4))
    Threads.@threads for vox in 1:tot_voxels
        if vec_mask[vox]
            x0[1] = maximum(vec_data[vox,:])
            x0[2] = vec_R2s_log[vox]
            x0[3] = vec_R2_log[vox]
            x0[4] = 1
            # model(xs,p) = sage_ns(echos,x0)
            fit = nlsfit(sage_ns_temp,echos,vec_data[vox,:],x0)
            # fit = Optim.minimizer(optimize(b -> sqerrorSAGE(b, echos,  DATA[IND[1],:,1]), x0))
            tmpOUT[vox,:] = fit;
            # vec_R2fit[vox,:] = fit[3];
        end
    end









    # x0 = Vector{Float64}(zeros(4))
    # @time Threads.@threads for dynamics in 1:nt
    #     for vox in 1:tot_voxels
    #         if vec_mask[vox]
    #             x0[1] = maximum(vec_data[vox,:,dynamics])
    #             x0[2] = vec_R2s_log[vox,dynamics]
    #             x0[3] = vec_R2_log[vox,dynamics]
    #             x0[4] = 1
    #             # model(xs,p) = sage_ns(echos,x0)
    #             fit = nlsfit(sage_ns,echos,vec_data[vox,:,dynamics],x0)
    #             # fit = Optim.minimizer(optimize(b -> sqerrorSAGE(b, echos,  DATA[IND[1],:,1]), x0))
    #             vec_R2sfit[vox,:] = fit[2];
    #             vec_R2fit[vox,:] = fit[3];
    #         end
    #     end
    #     println("Done with Fit $dynamics of $nt")
    # end