einstein_A_m1(k, L, S, J, I, F) = uconvert(u"s^-1", 64π^4 / 3h * k^3 * μ_0 / 4π * μ_B^2 * Float64(transition_matrix_element_m1_total(L, S, J, I, F)^2)) / (2F[1] + 1)

"iodine 127 magnetic diople transition line (;k, A, ν, λ, e)"
function line_I127(F1, F2)
    k = k_I127(F1, F2)
    A = A_I127(F1, F2)
    ν = uconvert(u"MHz", k * c_0)
    λ = uconvert(u"nm", 1 / k)
    e = uconvert(u"J", h * ν)
    (;k, A, ν, λ, e)
end

function k_I127(F1, F2)
    atom_name = "I127"
    atom_state1, atom_state2 = "e1", "g"
    df1 = zeeman_struc(atom_name, atom_state1)
	df2 = zeeman_struc(atom_name, atom_state2)
	k0 = df1[df1.F .== F1, :E][1] - df2[df2.F .== F2, :E][1]
end

function A_I127(F1, F2)
    L = (1, 1)
	S = (1/2, 1/2)
	J = (1/2, 3/2)
	I = (5/2, 5/2)
    k0 = k_I127(F1, F2)
    A0 = einstein_A_m1(k0, L, S, J, I, (F1, F2))
end

function σ0_I127(F1, F2, T, P, γ, M = ATOM_DATA[ATOM_DATA.Name .== "I127", :M][])
    A = A_I127(F1, F2)
    k0 = k_I127(F1, F2)
    ν0 = k0 * c_0
    νd = fwhm_doppler(ν0, M, T)
    νc = 2 * γ * P
    ν0, νd, νc = promote((ν0, νd, νc)...)
    σ0 = uconvert(u"cm^2*s^-1", A / k0^2 / 8π)
    voigt_peak = uconvert(u"s", profile_voigt(ν0; ν0 = ν0, νd = νd, νc = νc))
    σ0_peak = uconvert(u"cm^2", σ0 * voigt_peak)
    return (σ0, σ0_peak)
end