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

include("/mnt/c/Users/nicks/Documents/Github/The_MRI_toolbox/Dortch_additions_ToBe_deleted/models.jl")

# Load data from MATLAB
mat = matread("/mnt/c/Users/nicks/Documents/Github/The_MRI_toolbox/Dortch_additions_ToBe_deleted/matData.mat")
xmat = mat["x"]
p0mat = vec(mat["p0"])
ynmat = mat["yn"]
kmfmat = mat["kmf"]
Smmat = mat["Sm"]

println(size(ynmat))

# Define model - last argument is required for fixed kmf value, other optional kwargs can be defined
# Also we should be able to create another method that does not supply kmf for full model fitting

model(x,p) = SIR_Mz0(x,p,kmfmat,Sm=Smmat,R1m=NaN,mag=true) # use this to make sure we are using same values as MATLAB

function f()
    begin
        pj2 = zeros(size(ynmat))
        Threads.@threads for k in 1:size(ynmat,2)
            fit = curve_fit(model, xmat, ynmat[:,k], p0mat; autodiff=:finiteforward)
            pj2[:,k] = fit.param
        end
    end 
end


@time f() #3.366869 seconds (13.85 M allocations: 690.437 MiB, 6.12% gc time) 
@time f() #0.047446 seconds (393.96 k allocations: 32.638 MiB, 15.25% gc time)

# # Save to matlab for comparison
# matwrite("JuliaData.mat", Dict(
#                "pj1" => pj1,
#                "pj2" => pj2
#        ); compress = true)
