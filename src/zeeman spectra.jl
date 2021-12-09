"""

    zeeman_spec(atom, l1, l2, kx, T, νp, BF)

Compute the hyperfine spectra of several alkali atoms and iodine atom under Zeeman effect.
"""
function zeeman_spec(atom::Atom,
    l1::Int, l2::Int,
    kx::Vector{<:Wavenumber},
    T::Unitful.AbsoluteScaleTemperature, νp::Frequency,
    term::String,
    BF::BField = 0.0u"Gauss")

    level1 = l1 == 0 ? atom.ground : atom.excited[l1]
    (; L, S, J) = level1.state
    L1, S1, J1 = L, S, J
    level2 = l2 == 0 ? atom.ground : atom.excited[l2]
    (; L, S, J) = level2.state
    L2, S2, J2 = L, S, J
    df1 = zeeman_struc(atom, l1, BF)
    df2 = zeeman_struc(atom, l2, BF)
    n1 = nrow(df1)
    n2 = nrow(df2)
    df = DataFrame(
        F1 = HalfInteger[],
        MF1 = HalfInteger[],
        F2 = HalfInteger[],
        MF2 = HalfInteger[],
        q = HalfInteger[],
        k0 = Wavenumber[],
        c0 = Float64[],
        c = Vector{Float64}[],
        a = Frequency[],
        σ0 = Area[],
        σ = Vector{Area}[]
    )
    order = parse(Int, term[end])
    νx = uconvert.(unit(νp), kx * 𝑐)
    for i = 1:n1, j = 1:n2
        F1 = df1[i, :F]
        F2 = df2[j, :F]
        MF1 = df1[i, :MF]
        MF2 = df2[j, :MF]
        q = MF1 - MF2
        k0 = abs(df1[i, :E] - df2[j, :E])
        abs(q) > order && continue
        ## radiation matrix element
        ket1 = df1[i, :Ket2]
        ket2 = df2[j, :Ket2]
        c0 = relative_transitionME(ket1', ket2, order)^2
        ## line profile
        ν0 = uconvert(unit(νp), k0 * 𝑐)
        νd = fwhm_doppler(ν0, atom.M, T)
        peak = profile_voigt(ν0; ν0 = ν0, νd = νd, νp = νp)
        lineshape = profile_voigt.(νx; ν0 = ν0, νd = νd, νp = νp)
        c = c0 * ustrip.(lineshape)
        ## Einstein coefficient between two hfs state
        if term == "E1"
            ME = reducedME_E1(L1, S1, J1, L2, S2, J2)
        elseif term == "M1"
            ME = reducedME_M1(L1, S1, J1, L2, S2, J2)
        end
        a = aᵢⱼ(k0, c0 * ME^2)
        ## transition cross section
        σ = σᵢⱼ.(k0, a, lineshape)
        σ0 = σᵢⱼ(k0, a, peak)
        ## df
        push!(df, (F1, MF1, F2, MF2, q, k0, c0, c, a, σ0, σ))
    end
    sort!(df, [:F1, :F2], rev = true)
    return df
end