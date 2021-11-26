"""

    zeeman_struc(atom, level = 0, BF = 0u"Gauss")

Compute the hyperfine structure of several alkali atoms and iodine atom under Zeeman effect.
"""
function zeeman_struc(atom::Atom, level_num::Int = 0, BF::BField = 0u"Gauss")
    level = level_num == 0 ? atom.ground : atom.excited[level_num]
    (;L, S, J) = level.state
    I = atom.I
    (;A, B) = level.hfc
    basis_df = basis_hfs(L, S, J, I, couple = false)
    T_Ket1 = KetVec{Float64, UncoupledHyperfineStructureState{L, S, J, I}}
    T_Ket2 = KetVec{Float64, HyperfineStructureState{L, S, J, I}}
    df = DataFrame(
                    F = HalfInt[],
                    MF = HalfInt[],
                    ES = Wavenumber[],
                    E = Wavenumber[],
                    Ket1 = T_Ket1[],
                    Ket2 = T_Ket2[],
                )
    for MF in basis_df.MF
        basis1 = basis_df[basis_df.MF .== MF, :basis1][]
        basis2 = basis_df[basis_df.MF .== MF, :basis2][]
        gl = 1 - 𝑚e / atom.M
        h = hamiltonian_total(basis1, A, B, BF, μB = 𝜇B, μN = 𝜇N, gl = gl, gs = gS, gI = atom.gI)
        vals, vecs1 = diagonal(h)
        if A >= zero(A)&&B >= zero(B)
            reverse!(vals)
            reverse!(vecs1)
        end
        vecs2 = [basistransform(v, basis2) for v in vecs1]
        for i in 1:length(vals)
            vec1 = vecs1[i]
            vec2 = vecs2[i]
            val = vals[i]
            F = J + I - i + 1
            # a, p = findmax(vec2.c .^ 2)
            # F = vec2[p].s.F
            # a <= 0.5 && @warn "BF is so strong that hyperfine structure is substantially destructed." a F MF
            push!(df, [F, MF, val, val+level.E, vec1, vec2])
        end
    end
    sort!(df, [:F, :MF], rev=true)
    return df
end