function k_I127(F1, F2, MF1, MF2, BF = 0u"Gauss")
    df1 = zeeman_struc(I127, 0, BF)
	df2 = zeeman_struc(I127, 1, BF)
    E1 = filter(row -> row.F == F1 && row.MF == MF1, df1).E[]
    E2 = filter(row -> row.F == F2 && row.MF == MF2, df2).E[]
	return E2 - E1
end
k_I127(F1, F2) = k_I127(F1, F2, F1, F2, 0u"Gauss")

function A_I127(F1, F2)
    L = (1, 1)
	S = (1/2, 1/2)
	J = (3/2, 1/2)
	I = (5/2, 5/2)
    k0 = k_I127(F1, F2)
    return einstein_A_M1(k0, L, S, J, I, (F1, F2))
end

function σ0_I127(F1, F2; T::Unitful.AbsoluteScaleTemperature, P::Unitful.Pressure, γ = 5u"MHz/Torr")
    M = I127.M
    A = A_I127(F1, F2)
    k0 = k_I127(F1, F2)
    ν0 = k0 * 𝑐
    νd = fwhm_doppler(ν0, M, T)
    νp = 2 * γ * P
    σ0 = A / k0^2 / 8π |> u"cm^2*s^-1"
    voigt_peak = profile_voigt(ν0; ν0 = ν0, νd = νd, νp = νp) |> u"s"
    σ0_peak = σ0 * voigt_peak |> u"cm^2"
    return σ0_peak
end

function line_I127(F1, F2)
    k = k_I127(F1, F2)
    A = A_I127(F1, F2)
    ν = k * 𝑐 |> u"MHz"
    λ = 1 / k |> u"nm"
    e = ℎ * ν |> u"J"
    (;k, A, ν, λ, e)
end