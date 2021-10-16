"""
    zeeman_spec(atom_name::String, atom_state1::String, atom_state2::String, 
                kx::Vector{<:Unitful.Wavenumber},
                Temperature::Unitful.Temperature, vc::Unitful.Frequency,
                BF::Unitful.BField = 0u"Gauss")

Compute the hyperfine spectra of several alkali atoms and iodine atom under Zeeman effect.

# Arguments
- `atom_name::String`: "Li6", "K39", "K41", "Rb85", "Rb87", "Cs133", "I127"
- `atom_state::String`: "g","e1","e2" for alkali, "g","e1" for "I127"
- `kx::Vector{Unitful.Wavenumber`: spectra range
- `Temperature::Unitful.Temperature`: gas temperature, related to doppler broadening
- `vc::Unitful.Frequency`: gas pressure, related to collisional broadening
- `BF::Unitful.BField`: magnetic field strength

# Outputs
- `df_E::DataFrame`: F MF ES E
- `df_CS::DataFrame`: F1 MF1 F2 MF2 CrossSection-P CrossSection-S

For alkali atoms ,
- ground state : 2S1/2
- excited state 1 : 2P1/2(D1 line)
- excited state 2 : 2P3/2(D2 line)
For 127I,
- ground state : 2P3/2
- excited state 1 : 2P1/2
"""
function zeeman_spec(atom_name::String, atom_state1::String, atom_state2::String, 
                     kx::AbstractVector{<:Unitful.Wavenumber},
                     Temperature::Unitful.Temperature, νp::Unitful.Frequency,
                     BF::Unitful.BField = 0.0u"Gauss")
    atom_name in ATOM || error("zeeman_spec: atom_name=$(atom_name)")
    if atom_name in IODINE
        zeeman_spec_iodine(atom_name, atom_state1, atom_state2, kx, Temperature, νp, BF)
    elseif atom_name in ALKALI
        zeeman_spec_alkali(atom_name, atom_state1, atom_state2, kx, Temperature, νp, BF)
    end
end

function zeeman_spec_iodine(atom_name, atom_state1, atom_state2, kx, Temperature, νp, BF)
    if Set([atom_state1, atom_state2]) != Set(["g", "e1"])
        error("$atom_state1 -> $atom_state2 is not defined for $atom_name.")
    end
    #   1. compute the energy level split and quantum state functions
    df1_E::DataFrame, df1_coup::DataFrame, I1::HalfInt, J1::HalfInt, L1::HalfInt, S1::HalfInt, M1::typeof(1.0u"u") = 
        zeeman_struc(atom_name, atom_state1, BF; option = 2)
    df2_E::DataFrame, df2_coup::DataFrame, I2::HalfInt, J2::HalfInt, L2::HalfInt, S2::HalfInt, M2::typeof(1.0u"u") = 
        zeeman_struc(atom_name, atom_state2, BF; option = 2)

    #   2. compute the relative cross section between each state
    level1 = nrow(df1_E)
    level2 = nrow(df2_E)
    levelt = level1 * level2
    df_CS = DataFrame(
                F1 = Vector{HalfInteger}(undef, levelt), 
                MF1 = Vector{HalfInteger}(undef, levelt),
                F2 = Vector{HalfInteger}(undef, levelt),
                MF2 = Vector{HalfInteger}(undef, levelt),
                polarization = Vector{String}(undef, levelt),
                k0 = Vector{Quantity}(undef, levelt),
                a = Vector{Quantity}(undef, levelt),
                σ0 = Vector{Quantity}(undef, levelt),
                σ = Vector{Vector{Quantity}}(undef, levelt)
                )
    Threads.@threads for i in 1 : nrow(df1_E)
        for j in 1 : nrow(df2_E)
            k::typeof(1.0u"cm^-1") = abs(df1_E[i, :E] - df2_E[j, :E])
            MF1::HalfInt = df1_E[i, :MF]
            MF2::HalfInt = df2_E[j, :MF]
            q = MF1 - MF2
            if q == 0
                polarization = "P"
            elseif q == -1 || q == 1
                polarization = "S"
            else
                polarization = "M1 forbidden"
            end
            sum = 0.0
            if q in HalfInt.((0, 1, -1))
                for c1 in 3 : size(df1_coup, 2)
                    for c2 in 3 : size(df2_coup, 2)
                        if df1_coup[i, c1]*df2_coup[j, c2] > eps()
                            F1c = I1 + J1 + 3 - c1
                            F2c = I2 + J2 + 3 - c2
                            L = (L1, L2)
                            S = (S1, S2)
                            J = (J1, J2)
                            I = (I1, I2)
                            F = (F1c, F2c)
                            MF = (MF1, MF2)
                            sum += df1_coup[i, c1] * df2_coup[j, c2] * 
                                transition_matrix_element_m1(L, S, J, I, F, MF, q)
                        end
                    end
                end
            end
            F1 = df1_E[i, :F]::HalfInt
            F2 = df2_E[j, :F]::HalfInt
            ac::typeof(1.0u"s^-1") = uconvert(u"s^-1", 64π^4/3/h*k^3*μ_0/4π*μ_B^2)
            a = sum^2 < eps() ? zero(ac) : ac * sum^2
            ν0 = k * c_0
            νd = fwhm_doppler(ν0, M1, Temperature)
            σ0 = uconvert.(u"cm^2 * s^-1", a / 8π / k^2)
            ν0, νd, νp = promote((ν0, νd, νp)...)
            νx = uconvert.(unit(ν0), kx * c_0)
            σ = σ0 * profile_voigt.(νx; ν0 = ν0, νd = νd, νp = νp)
            σ = uconvert.(u"cm^2", σ)
            n = (i - 1) * level2 + j
            df_CS[n, :] = F1, MF1, F2, MF2, polarization, k, a, σ0, σ
        end
    end
    df_E = append!(df1_E, df2_E)
    return df_E, df_CS
end

function zeeman_spec_alkali(atom_name, atom_state1, atom_state2, kx, Temperature, νp, BF)
    if Set([atom_state1, atom_state2]) != Set(["g", "e1"])&&
            Set([atom_state1, atom_state2]) != Set(["g", "e2"])
        error("$atom_state1 -> $atom_state2 is not defined for $atom_name.")
    end
    #   1. compute the energy level split and quantum state functions
    df1_E::DataFrame, df1_coup::DataFrame, I1::HalfInt, J1::HalfInt, L1::HalfInt, S1::HalfInt, M1::typeof(1.0u"u") = 
        zeeman_struc(atom_name, atom_state1, BF; option = 2)
    df2_E::DataFrame, df2_coup::DataFrame, I2::HalfInt, J2::HalfInt, L2::HalfInt, S2::HalfInt, M2::typeof(1.0u"u") = 
        zeeman_struc(atom_name, atom_state2, BF; option = 2)
    #   2. compute the relative cross section between each state
    ar0 = 0
    for F1 = J1+I1:-1:abs(J1-I1)
        for F2 = J2+I2:-1:abs(J2-I2)
            cJIF = uncoup_T1k(J1, I1, F1, J2, I2, F2, 1)^2
            cLSJ = uncoup_T1k(L1, S1, J1, L2, S2, J2, 1)^2
            ar0 = cJIF * cLSJ > ar0 ? cJIF * cLSJ : ar0
        end
    end
    level1 = nrow(df1_E)
    level2 = nrow(df2_E)
    levelt = level1 * level2
    df_CS = DataFrame(
                F1 = Vector{HalfInteger}(undef, levelt), 
                MF1 = Vector{HalfInteger}(undef, levelt),
                F2 = Vector{HalfInteger}(undef, levelt),
                MF2 = Vector{HalfInteger}(undef, levelt),
                polarization = Vector{String}(undef, levelt),
                k0 = Vector{Quantity}(undef, levelt),
                ar = Vector{Float64}(undef, levelt),
                σr = Vector{Vector{Float64}}(undef, levelt)
                )
    Threads.@threads for i in 1 : level1
        for j in 1 : level2
            k::typeof(1.0u"cm^-1") = abs(df1_E[i, :E] - df2_E[j, :E])
            MF1::HalfInt = df1_E[i, :MF]
            MF2::HalfInt = df2_E[j, :MF]
            q = MF1 - MF2
            if q == 0
                polarization = "P"
            elseif q == -1 || q == 1
                polarization = "S"
            else
                polarization = "E1 forbidden"
            end
            sum = 0.0
            if q in HalfInt.((0, 1, -1))
                for c1 in 3 : size(df1_coup, 2)
                    for c2 in 3 : size(df2_coup, 2)
                        if df1_coup[i, c1]*df2_coup[j, c2] > eps()
                            F1c = I1 + J1 + 3 - c1
                            F2c = I2 + J2 + 3 - c2
                            L = (L1, L2)
                            S = (S1, S2)
                            J = (J1, J2)
                            I = (I1, I2)
                            F = (F1c, F2c)
                            MF = (MF1, MF2)
                            sum += df1_coup[i, c1] * df2_coup[j, c2] * 
                                transition_matrix_element_e1(L, S, J, I, F, MF, q)
                        end
                    end
                end
            end
            νd = fwhm_doppler(k*c_0, M1, Temperature)
            ar = sum^2 * 3 / ar0
            ν0 = k * c_0
            νx = kx * c_0
            ν0, νd, νp = promote((ν0, νd, νp)...)
            νx = uconvert.(unit(ν0), kx * c_0)
            σr = ar * ustrip.(NoDims, profile_voigt.(νx; ν0 = ν0, νd = νd, νp = νp)/
                            profile_voigt(ν0; ν0 = ν0, νd = νd, νp = νp))
            F1 = df1_E[i, :F]::HalfInt
            F2 = df2_E[j, :F]::HalfInt
            n = (i - 1) * level2 + j
            df_CS[n, :] = F1, MF1, F2, MF2, polarization, k, ar, σr
        end
    end
    df_E = append!(df1_E, df2_E)
    return df_E, df_CS
end