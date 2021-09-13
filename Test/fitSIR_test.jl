using Pkg	
try
    println("If is first time you ran the code. It will take a minute to precompile.")
    @eval using Statistics; 
    @eval using LsqFit;
    @eval using MAT
    Pkg.precompile()
catch e
    # not found; install and try loading again
    Pkg.add("Statistics")
    Pkg.add("LsqFit")
    Pkg.add("MAT")
    @eval using Statistics; 
    @eval using LsqFit;
    @eval using MAT
end

using NIfTI; 
using LsqFit;
using Printf
using MAT

include("/mnt/c/Users/nicks/Documents/Github/The_SIR-qMT_toolbox/Test/models.jl")
include("/mnt/c/Users/nicks/Documents/Github/The_SIR-qMT_toolbox/utils.jl")

# Load data from MATLAB
mat = matread("/mnt/c/Users/nicks/Downloads/matData.mat")
xmat = mat["x"]
p0mat = vec(mat["p0"])
ynmat = mat["yn"]
kmfmat = mat["kmf"]
Smmat = mat["Sm"]

println(size(ynmat))

# Define model - last argument is required for fixed kmf value, other optional kwargs can be defined
# Also we should be able to create another method that does not supply kmf for full model fitting

ti = xmat[:,1]
td = xmat[:,2]
# model(x,p) = SIR_Mz0(x,p,kmfmat,Sm=Smmat,R1m=NaN,mag=true) # use this to make sure we are using same values as MATLAB

function SIR_Mz0_v2(ti::Vector{Float64}, p::Vector{Float64},td::Vector{Float64}, kmf::Float64; 
    Sm::Float64=0.83, R1m::Float64=NaN, mag::Bool=true)

    # Extract ti and td values from x
    # ti = x[:,1]
    # td = x[:,2]

    # Define model parameters based on p
    pmf = p[1]
    R1f = p[2]
    Sf  = p[3]
    Mf∞ = p[4]

    # Define R1m based on user-defined value (=R1f when set to NaN)
    if isnan(R1m)
        R1m = R1f
    end

    # Define kfm based on kmf and pmf (assuming mass balance)
    kfm = kmf*pmf

    # Apparent rate constants
    ΔR1 = sqrt((R1f-R1m+kfm-kmf)^2.0 + 4.0*kfm*kmf)
    R1⁺ = (R1f + R1m + kfm + kmf + ΔR1) / 2.0
    R1⁻ = R1⁺ - ΔR1

    # Component amplitudes for td terms
    bf_td⁺ = -(R1f - R1⁻) / ΔR1
    bf_td⁻ =  (R1f - R1⁺) / ΔR1
    bm_td⁺ = -(R1m - R1⁻) / ΔR1
    bm_td⁻ =  (R1m - R1⁺) / ΔR1

    # Loop over ti/td values
    # make this a new function
    M = similar(ti)
    for k in 1:length(td)

        # Signal recovery during td
        E_td⁺ = exp(-R1⁺*td[k])
        E_td⁻ = exp(-R1⁻*td[k])
        Mf_td = bf_td⁺*E_td⁺ + bf_td⁻*E_td⁻ + 1.0
        Mm_td = bm_td⁺*E_td⁺ + bm_td⁻*E_td⁻ + 1.0

        # Component amplitude terms for ti terms
        a = Sf*Mf_td - 1.0
        b = (Sf*Mf_td - Sm*Mm_td) * kfm
        bf_ti⁺ =  (a*(R1f-R1⁻) + b) / ΔR1
        bf_ti⁻ = -(a*(R1f-R1⁺) + b) / ΔR1

        # Signal recovery during ti
        M[k] = (bf_ti⁺*exp(-R1⁺*ti[k]) + bf_ti⁻*exp(-R1⁻*ti[k]) + 1.0) * Mf∞

        # Take the magnitude of the signal
        if mag
            M[k] = abs(M[k])
        end
    end

    # Return signal
    # return Vector{Float64}(M)
    return oftype(M,ti) # drammatic speed up
end
model(x,p) = SIR_Mz0_v2(ti,p,td,kmfmat,Sm=Smmat,R1m=NaN,mag=true) # use this to make sure we are using same values as MATLAB

function f()
    begin
        pj2 = zeros(size(ynmat))
        Threads.@threads for k in 1:size(ynmat,2)
            fit = curve_fit(model, ti, ynmat[:,k], p0mat; autodiff=:finiteforward)
            pj2[:,k] = fit.param
        end
    end 
end

@time f() #3.366869 seconds (13.85 M allocations: 690.437 MiB, 6.12% gc time) 
@time f() #0.047446 seconds (393.96 k allocations: 32.638 MiB, 15.25% gc time)

temp = niread("/mnt/c/Users/nicks/Documents/MRI_data/PING_brains/TestRetest/output_20210618/proc_20210618/brain.nii.gz");
temp.header
temp2=deepcopy(temp)
DATA = Array{Float32}(temp.raw);
d=voxel_size(temp.header)
voxset = (d[1],d[2],d[3])
temp2 = NIVolume(DATA; voxel_size=(voxset))
aff=getaffine(temp)
aff = [d[1] 0 0 -aff[1,end];
    0 d[2] 0 -aff[2,end];
    0 0 d[3] aff[3,end];
    0 0 0 1]  

setaffine(temp2.header,aff)
niwrite("/mnt/c/Users/nicks/Desktop/temp.nii.gz",temp2)


# # Save to matlab for comparison
# matwrite("JuliaData.mat", Dict(
#                "pj1" => pj1,
#                "pj2" => pj2
#        ); compress = true)
