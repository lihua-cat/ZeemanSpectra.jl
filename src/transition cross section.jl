function einstein_A_M1(k::Wavenumber,
                       L::NTuple{2},
                       S::NTuple{2},
                       J::NTuple{2},
                       I::NTuple{2}, 
                       F::NTuple{2})
    c1 = 64π^4 / 3ℎ * k^3 * 𝜇0 / 4π * μ_B^2
    c2 = uncoup_T1(J[1], I[1], F[1], J[2], I[2], F[2], 1) |> Float64
    ME = reducedME_M1(L[1], S[1], J[1], L[2], S[2], J[2])
    c3 = 2F[2] + 1
    return c1 * (c2 * ME)^2 / c3 |> u"s^-1"
end

function einstein_A_M1(k, s1, s2)
    L1, S1, J1, I1, F1 = s1
    L2, S2, J2, I2, F2 = s2
    L = (L1, L2)
    S = (S1, S2)
    J = (J1, J2)
    I = (I1, I2)
    F = (F1, F2)
    return einstein_A_M1(k, L, S, J, I, F)
end

function k_I127(F1, F2)
    df1 = zeeman_struc(I127, 0)
	df2 = zeeman_struc(I127, 1)
	return df2[df2.F .== F2, :E][1] - df1[df1.F .== F1, :E][1]
end

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