function k_I127(F1, F2, MF1, MF2, BF = 0u"Gauss")
    df1 = zeeman_struc(I127, 0, BF)
	df2 = zeeman_struc(I127, 1, BF)
    E1 = filter(row -> row.F == F1 && row.MF == MF1, df1).E[]
    E2 = filter(row -> row.F == F2 && row.MF == MF2, df2).E[]
	return E2 - E1
end
k_I127(F1, F2) = k_I127(F1, F2, F1, F2, 0u"Gauss")

function A_I127(F1, F2)
    L = 1
	S = 1/2
	J = (3/2, 1/2)
	I = 5/2
    k0 = k_I127(F1, F2)
    return einsteinA(k0, L, S, J[1], I, F1, L, S, J[2], I, F2, "M1")
end

"g = σ * (Nu - gl/gu * Nl)"
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

"g = σ * (Nᵢ - Nⱼ)"
function σm_I127(F1, F2, BF, p; T, P, γ)
	kx = collect(7602.2:0.0001:7603.8)u"cm^-1"
	df = zeeman_spec(I127, 0, 1, kx, T, 2P*γ, "M1", BF)
	dff = filter(row->row.F1==F1&&row.F2==F2, df)
	if p == "P"
		dfff = filter(row->row.q==0, dff)
		σ = sum(dfff.σ) * 3
	elseif p == "S"
		dfff = filter(row->abs(row.q)==1, dff)
		σ = sum(dfff.σ) / 2 * 3
	end
    σm = maximum(σ)
	return σm
end

function σr_ltp(F1, F2, B_sample, p; T, P, γ)
	σ0 = σ0_I127(F1, F2, T=T, P=P, γ=γ) * (2F2 + 1)
	σm = σm_I127.(F1, F2, B_sample*u"Gauss", p, T=T, P=P, γ=γ) / σ0
	ltp = LinearInterpolation(B_sample, σm)
	return ltp
end