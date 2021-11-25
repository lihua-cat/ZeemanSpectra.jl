function Hamiltonian(J1::T, I1::T, MJ1::T, MI1::T, 
                     J2::T, I2::T, MJ2::T, MI2::T, 
                     gJ::Float64, gI::Float64, 
                     A::Float64, B::Float64, 
                     BF::Number) where {T<:HalfInteger}
    (BF >= 0 && J1 == J2 && I1 == I2) || error("Hamiltonian: ")
    quantum_number = (J1, I1, MJ1, MI1, J2, I2, MJ2, MI2)
    Jz = angular(quantum_number..., ["Jz"])
    Iz = angular(quantum_number..., ["Iz"])
    IpJm = angular(quantum_number..., ["Ip", "Jm"])
    ImJp = angular(quantum_number..., ["Im", "Jp"])
    JzJpIzIm = angular(quantum_number..., ["Jz", "Jp", "Iz", "Im"])
    JzJpImIz = angular(quantum_number..., ["Jz", "Jp", "Im", "Iz"])
    JpJzIzIm = angular(quantum_number..., ["Jp", "Jz", "Iz", "Im"])
    JpJzImIz = angular(quantum_number..., ["Jp", "Jz", "Im", "Iz"])
    JzJmIzIp = angular(quantum_number..., ["Jz", "Jm", "Iz", "Ip"])
    JzJmIpIz = angular(quantum_number..., ["Jz", "Jm", "Ip", "Iz"])
    JmJzIzIp = angular(quantum_number..., ["Jm", "Jz", "Iz", "Ip"])
    JmJzIpIz = angular(quantum_number..., ["Jm", "Jz", "Ip", "Iz"])
    Ip2Jm2 = angular(quantum_number..., ["Ip", "Ip", "Jm", "Jm"])
    Im2Jp2 = angular(quantum_number..., ["Im", "Im", "Jp", "Jp"])
    c1 = (MJ1, MI1) == (MJ2, MI2) ? (3MJ2^2 - aa(J2)) * (3MI2^2 - aa(I2)) / 2 : 0
    Hzm = μB * BF * gJ * Jz - μN * BF * gI * Iz
    Hhf = A * (Iz * Jz + (IpJm + ImJp) / 2) +
          B / 2I2 / (2I2 - 1) / J2 / (2J2 - 1) * 
          (
              c1 + 3 / 4 * 
              (
                  JzJpIzIm + JzJpImIz + JpJzIzIm + JpJzImIz +
                  JzJmIzIp + JzJmIpIz + JmJzIzIp + JmJzIpIz +
                  Ip2Jm2 + Im2Jp2
              )
          )
    return Hzm + Hhf
end

"""
zeeman_struc(atom_name, atom_state, BF = 0u"Gauss")

Compute the hyperfine structure of several alkali atoms and iodine atom under Zeeman effect.

# Arguments
- `atom_name::String`: "Li6", "K39", "K41", "Rb85", "Rb87", "Cs133", "I127"
- `atom_state::String`: "g","e1","e2" for alkali, "g","e1" for "I127"
- `BF::Unitful.BField`: magnetic field strength

# Outputs
- `df::DataFrame`: F MF ES E |F,MF> |F-1,MF> ... |MJ,MI> ...

For alkali atoms,
- ground state : 2S1/2
- excited state 1 : 2P1/2(D1 line)
- excited state 2 : 2P3/2(D2 line)
For 127I,
- ground state : 2P3/2
- excited state 1 : 2P1/2
"""
function zeeman_struc(atom_name::String,
            atom_state::String,
            BF::Unitful.BField = 0u"Gauss";
            option::Int64 = 1)

    atom_name in ATOM || error("zeeman_struc: atom_name=$(atom_name)")
    if atom_name in ALKALI
        atom_state in ("g", "e1", "e2") || error("zeeman_struc: atom_state=$(atom_state)")
    elseif atom_name in IODINE
        atom_state in ("g", "e1") || error("zeeman_struc: atom_state=$(atom_state)")
    end
    BF = ustrip(u"Gauss", BF)
    #   1. Load data and physical constants
    row = findfirst(ATOM_DATA.Name .== atom_name)
    dfrow = ATOM_DATA[row, :]
    #   atomic weight
    Mu = dfrow.M
    M = ustrip(u"u", Mu)
    #   g-factor
    gI = dfrow.gI
    gL = 1 - me / M
    gs = 1.00115965218085 * 2 #https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.030801
    #   angular quantum number
    S = 1//2
    I = dfrow.I
    if atom_name in ALKALI
    if atom_state == "g"
    L = 0
    J = 1//2
    else
    L = 1
    J = atom_state == "e1" ? 1//2 : 3//2
    end
    elseif atom_name in IODINE
    L = 1
    J = atom_state == "g" ? 3//2 : 1//2
    end
    S, L, J, I = HalfInt.((S, L ,J, I))
    gJ = gL * (aa(J) - aa(S) + aa(L))/2aa(J) + 
        gs * (aa(J) + aa(S) - aa(L))/2aa(J)
    #   hyperfine interaction coefficients A, B
    A = ustrip(u"cm^-1", dfrow[Symbol("A"*atom_state)])
    B = ustrip(u"cm^-1", dfrow[Symbol("B"*atom_state)])
    #   reference energy level
    E0u = atom_state == "g" ? 0u"cm^-1" : 
    atom_state == "e1" ? dfrow[Symbol("D1")] : dfrow[Symbol("D2")]
    E0 = ustrip(u"cm^-1", E0u)

    #   2. Initial level and coe dataframe
    jx = min(I ,J)
    jd = max(I ,J)
    n_level = (2J + 1) * (2I + 1) |> Int # number of hyperfine splits
    n_state = 2jx + 1 |> Int   # max number of coupling states, |F, MF> -> n |MJ, MI>
    F_list = I + J : -1 : abs(I - J)
    MF_list = I + J : -1 : -(I + J)
    Mx_list = jx : -1 : -jx
    Vec_F = Vector{HalfInteger}(undef, n_level)
    Vec_MF = Vector{HalfInteger}(undef, n_level)
    Mat_uncoup = zeros(Float64, n_level, n_state)
    Mat_coup = zeros(Float64, n_level, n_state)
    n = 1
    for F in F_list
    MF_vec = collect(F:-1:-F) 
    Vec_F[n : n + Int(2F + 1) - 1] .= F
    Vec_MF[n : n + Int(2F + 1) - 1] .= MF_vec
    n += Int(2F + 1)
    end

    #   3. degenerate perturbation theory: solve det(H'-eigE) = 0
    #   special case
    if B == zero(B) && J == 1//2
        x = 2(gJ * μB + gI * μN) * BF / (2I + 1)A
        sq = @. √(x^2 + 4x * Vec_MF / (2I + 1) + 1)
        sq[Vec_MF .== -I - 1/2] .= 1 - x
        Vec_ES = @. -A/4 - gI * μN * BF * Vec_MF + 2(Vec_F - I) * (2I + 1)A / 4 * sq
        dd = @. 1/2 * A * √((I - Vec_MF + 1/2) * (I + 1/2 + Vec_MF))
        ddd = @. Vec_ES - 1/2 * A * (Vec_MF - 1/2)- 1/2 * gJ * μB * BF + 
        (Vec_MF - 1/2) * gI * μN * BF
        d = ddd[abs.(Vec_MF) .!= I+J] ./ dd[abs.(Vec_MF) .!= I+J]
        Mat_uncoup[abs.(Vec_MF) .!= I+J, :] .= [one.(d) d] ./ sqrt.(1 .+ d.^2)
        Mat_uncoup[Vec_MF .== I+J, :] .= [1 0]
        Mat_uncoup[Vec_MF .== -(I+J), :] .= [0 1]
        Vec_E = Vec_ES .+ E0
    else
        #   general case
        Vec_ES = Vector{Float64}(undef, n_level)
        for MF in MF_list
            Md_list = Vector(filter(x -> abs(x) <= jd, MF .- Mx_list))
            Mxx_list = MF .- Md_list
            Hmn = Matrix{Float64}(undef, length(Md_list), length(Md_list))
            Threads.@threads for i in CartesianIndices(Hmn)
                    if J >= I
                        MJ1 = Md_list[i[1]]
                        MI1 = Mxx_list[i[1]]
                        MJ2 = Md_list[i[2]]
                        MI2 = Mxx_list[i[2]]
                    else
                        MJ1 = Mxx_list[i[1]]
                        MI1 = Md_list[i[1]]
                        MJ2 = Mxx_list[i[2]]
                        MI2 = Md_list[i[2]]
                    end
                Hmn[i] = Hamiltonian(J, I, MJ1, MI1, J, I, MJ2, MI2, gJ, gI, A, B, BF)
            end
            vals, vecs = eigen(Hmn, sortby = x -> -x)
            vecs = transpose(vecs)
            vecs = sign.(vecs[:, 1]) .* vecs
            Vec_ES[Vec_MF .== MF] = vals
            n1 = Int(1 + jx - Mxx_list[1])
            n2 = Int(n1 + length(Mxx_list) - 1)
            Mat_uncoup[Vec_MF .== MF, n1:n2] = vecs
        end
        Vec_E = Vec_ES .+ E0
    end
    df_E = DataFrame(F = Vec_F, MF = Vec_MF, 
                     ES = Vec_ES*1u"cm^-1", E = Vec_E*1u"cm^-1")
    df_uncoup = DataFrame(F = Vec_F, MF = Vec_MF)
    s1 = J >= I ? "MI" : "MJ"
    for i in 1:length(Mx_list)
        s2 = string(Mx_list[i])
        s = s1 * "=" * s2
        df_uncoup.B = Mat_uncoup[:, i]
        rename!(df_uncoup, :B => s)
    end
    #   4. basis transformation
    Threads.@threads for i in 1 : n_level
        for j in 1 : n_state
            F = I + J + 1 - j
            for k = 1 : n_state
                if abs(Mat_uncoup[i, k]) > eps()&&F >= abs(Vec_MF[i])
                    Mx = jx + 1 - k
                    @inbounds Mat_coup[i, j] += Mat_uncoup[i, k] * 
                        clebschgordan(jx, Mx, jd, Vec_MF[i] - Mx, F)
                end
            end
        end
    end
    Mat_coup[abs.(Mat_coup) .< eps()] .= 0.0
    df_coup = DataFrame(F = Vec_F, MF = Vec_MF)
    for i in 1 : length(F_list)
        s = "F=" * string(F_list[i])
        df_coup.B = Mat_coup[:, i]
        rename!(df_coup, :B => s)
    end
    #   5.  done
    if option == 1 
        return hcat(df_E, df_coup[!, 3:end], df_uncoup[!, 3:end])
    elseif option == 2
        return df_E, df_coup, I, J, L, S, Mu
    end
end