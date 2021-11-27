aᵢⱼ(k, s) = 64π^4 / 3ℎ * k^3 * s |> u"s^(-1)"

function einsteinA(k::Wavenumber,
                   L::NTuple{2},
                   S::NTuple{2},
                   J::NTuple{2},
                   I::NTuple{2},
                   F::NTuple{2},
                   term::String)
    c = uncoup_T1(J[1], I[1], F[1], J[2], I[2], F[2], 1) |> Float64
    if term == "E1"
        ME = reducedME_E1(L[1], S[1], J[1], L[2], S[2], J[2])
    elseif term == "M1"
        ME = reducedME_M1(L[1], S[1], J[1], L[2], S[2], J[2])
    end
    g = 2F[2] + 1
    s = (c * ME)^2
    gA = aᵢⱼ(k, s)
    A = gA / g |> u"s^-1"
    return A
end

function einsteinA(k, s1, s2, term)
    L1, S1, J1, I1, F1 = s1
    L2, S2, J2, I2, F2 = s2
    L = (L1, L2)
    S = (S1, S2)
    J = (J1, J2)
    I = (I1, I2)
    F = (F1, F2)
    return einsteinA(k, L, S, J, I, F, term)
end

