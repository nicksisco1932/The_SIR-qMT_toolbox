
# using NIfTI; 
# using LsqFit;
# using Printf
using ArgParse;

# include("./utils.jl")
function commandline()
        

    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "SIR_nii_1"
        required = true
        "SIR_nii_2"
        required = true
        "SIR_nii_3"
        required = true
        "SIR_nii_4"
        required = true
        "SIR_nii_brainMask"
        required = true
    end

    println(parse_args(settings))

    # parsed_args = SIR_parse_commandline();
    for (out, val) in parse_args(settings)
        println(" $out => $val")
    end
    return parse_args(settings)
end

a = commandline()

base = a["SIR_nii_1"]

