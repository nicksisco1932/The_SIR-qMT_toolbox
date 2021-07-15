function SIR_Mz0(x::Matrix{Float64},p::Vector{Float64}, kmf::Float64;
    Sm::Float64=0.83, R1m::Float64=NaN, mag::Bool=true)

    # Extract ti and td values from x
    ti = x[:,1]
    td = x[:,2]

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
    return M
end
