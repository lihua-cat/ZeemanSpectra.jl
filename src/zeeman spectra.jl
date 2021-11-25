"""

    zeeman_spec(atom, l1, l2, kx, T, νp, BF)

Compute the hyperfine spectra of several alkali atoms and iodine atom under Zeeman effect.
"""
function zeeman_spec(atom::Atom, 
                     l1::Int, l2::Int, 
                     kx::Vector{<:Wavenumber},
                     T::Unitful.AbsoluteScaleTemperature, νp::Frequency,
                     BF::BField = 0.0u"Gauss")
    df1 = zeeman_struc(atom, l1, BF)
    df2 = zeeman_struc(atom, l2, BF)
    n1 = nrow(df1)
    n2 = nrow(df2)
    nt = n1 * n2
    df = DataFrame(
                    F1 = Vector{HalfInteger}(undef, nt), 
                    MF1 = Vector{HalfInteger}(undef, nt),
                    F2 = Vector{HalfInteger}(undef, nt),
                    MF2 = Vector{HalfInteger}(undef, nt),
                    q = Vector{HalfInteger}(undef, nt),
                    k0 = Vector{Wavenumber}(undef, nt),
                    c0 = Vector{Float64}(undef, nt),
                    c = Vector{Vector{Float64}}(undef, nt),
                    a = Vector{Frequency}(undef, nt),
                    σ0 = Vector{Area}(undef, nt),
                    σ = Vector{Vector{Area}}(undef, nt)
                  )
    Threads.@threads for i in 1 : n1
        for j in 1 : n2
            F1 = df1[i, :F]
            F2 = df2[j, :F]
            MF1 = df1[i, :MF]
            MF2 = df2[j, :MF]
            q = MF1 - MF2
            k0 = abs(df1[i, :E] - df2[j, :E])
            ## radiation matrix element
            ket1 = df1[i, :Ket2]
            ket2 = df2[j, :Ket2]
            c0 = relative_transition_intensity(ket1', ket2)^2
            ## line profile
            ν0 = uconvert(unit(νp), k0 * 𝑐)
            νd = fwhm_doppler(ν0, atom.M, T)
            νx = uconvert.(unit(νp), kx * 𝑐)
            peak = profile_voigt(ν0; ν0 = ν0, νd = νd, νp = νp)
            lineshape = profile_voigt.(νx; ν0 = ν0, νd = νd, νp = νp)
            c = c0 * ustrip.(lineshape)
            ## Einstein coefficient between two hfs state
            C1 = 64π^4 / 3ℎ * k0^3 * 𝜇0 / 4π * μ_B^2
            a = C1 * c0 |> u"s^(-1)"
            ## transition cross section
            C2 = 1 / k0^2 / 8π
            σ = C2 * a * lineshape .|> u"cm^2"
            σ0 = C2 * a * peak |> u"cm^2"
            ## df
            n = (i - 1) * n2 + j
            df[n, :] = F1, MF1, F2, MF2, q, k0, c0, c, a, σ0, σ
        end
    end
    return df
end