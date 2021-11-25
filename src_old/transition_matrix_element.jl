@doc raw"""

    transition_matrix_element_m1(L, S, J, I, F, MF, q)

compute MD transition matrix element $\langle \alpha JIFM_F|\mathrm{(J+S)}_q|\alpha' J' I' F' M_F'\rangle$.
"""
function transition_matrix_element_m1(L::NTuple{2, T}, S::NTuple{2, T}, J::NTuple{2, T}, 
        I::NTuple{2, T}, F::NTuple{2, T}, MF::NTuple{2, T}, q::T) where {T <: HalfInteger}
    #   <FMF|T^1_q|F'MF'> = c1 * <F||T^1||F'>
    c1 = we_Tkq(F[1], F[2], MF[1], MF[2], 1, q)
    #   <JIF||T^1||J'I'F'> = c2 * <J||T^1||J>
    c2 = uncoup_T1k(J[1], I[1], F[1], J[2], I[2], F[2], 1)
    #   <LSJ||(J+S)||L'S'J'> = <LSJ||J||L'S'J'> + <LSJ||S||L'S'J'>
    rme_LSJ_J = rme_J(J[1], J[2])
    #   <LSJ||S||L'S'J'> = c3 * <S||S||S>
    c3 = uncoup_T2k(L[1], S[1], J[1], L[2], S[2], J[2], 1)
    rme_S = rme_J(S[1], S[2])
    return c1 * c2 * RationalRoot(rme_LSJ_J + c3 * rme_S)
end
function transition_matrix_element_m1(L::NTuple{2,}, S::NTuple{2,}, J::NTuple{2,}, 
        I::NTuple{2,}, F::NTuple{2,}, MF::NTuple{2,}, q::Real)
    # check quantum number isa halfintrger ?
    all(ishalfinteger, (L..., S..., J..., I..., F..., MF..., q...)) || error("transition_matrix_element_m1: ")
    L = HalfInt.(L)
    S = HalfInt.(S)
    J = HalfInt.(J)
    I = HalfInt.(I)
    F = HalfInt.(F)
    MF = HalfInt.(MF)
    q = HalfInt.(q)
    transition_matrix_element_m1(L, S, J, I, F, MF, q)
end


@doc raw"""

    transition_matrix_element_e1(L, S, J, I, F, MF, q)

compute `C` of ED transition matrix element $\langle \alpha JIFM_F|\mathrm{r}_q|\alpha' J' I' F' M_F'\rangle$.
```math
\langle \alpha JIFM_F|\mathrm{r}_q|\alpha' J' I' F' M_F'\rangle = C * \langle \alpha L\left\|\mathrm{r}\right\|\alpha L'\rangle
```
"""
function transition_matrix_element_e1(L::NTuple{2, T}, S::NTuple{2, T}, J::NTuple{2, T}, 
        I::NTuple{2, T}, F::NTuple{2, T}, MF::NTuple{2, T}, q::T) where {T <: HalfInteger}
    #   <FMF|T^1_q|F'MF'> = c1 * <F||T^1||F'>
    c1 = we_Tkq(F[1], F[2], MF[1], MF[2], 1, q)
    #   <JIF||T^1||J'I'F'> = c2 * <J||T^1||J>
    c2 = uncoup_T1k(J[1], I[1], F[1], J[2], I[2], F[2], 1)
    #   <LSJ||r||L'S'J'> = c3 * <L||r||L>
    #   <L||r||L'> depends on ψ(r), radical wave function.
    c3 = uncoup_T1k(L[1], S[1], J[1], L[2], S[2], J[2], 1)
    return c1 * c2 * c3
end
function transition_matrix_element_e1(L::NTuple{2,}, S::NTuple{2,}, J::NTuple{2,}, 
    I::NTuple{2,}, F::NTuple{2,}, MF::NTuple{2,}, q::Real)
    # check quantum number isa halfintrger ?
    all(ishalfinteger, (L..., S..., J..., I..., F..., MF..., q...)) || error("transition_matrix_element_e1: ")
    L = HalfInt.(L)
    S = HalfInt.(S)
    J = HalfInt.(J)
    I = HalfInt.(I)
    F = HalfInt.(F)
    MF = HalfInt.(MF)
    q = HalfInt.(q)
    transition_matrix_element_e1(L, S, J, I, F, MF, q)
end



@doc raw"""

    transition_matrix_element_m1_total(L, S, J, I, F)

compute MD transition matrix element $\langle \alpha JIF\left\|\mathrm{(J+S)}\right\|\alpha' J' I' F'\rangle$.
"""
function transition_matrix_element_m1_total(
        L::NTuple{2, T}, S::NTuple{2, T}, J::NTuple{2, T}, 
        I::NTuple{2, T}, F::NTuple{2, T}) where {T <: HalfInteger}
    #   <JIF||T^1||J'I'F'> = c2 * <J||T^1||J>
    c2 = uncoup_T1k(J[1], I[1], F[1], J[2], I[2], F[2], 1)
    #   <LSJ||(J+S)||L'S'J'> = <LSJ||J||L'S'J'> + <LSJ||S||L'S'J'>
    rme_LSJ_J = rme_J(J[1], J[2])
    #   <LSJ||S||L'S'J'> = c3 * <S||S||S>
    c3 = uncoup_T2k(L[1], S[1], J[1], L[2], S[2], J[2], 1)
    rme_S = rme_J(S[1], S[2])
    return c2 * RationalRoot(rme_LSJ_J + c3 * rme_S)
end
function transition_matrix_element_m1_total(L::NTuple{2,}, S::NTuple{2,}, J::NTuple{2,}, 
    I::NTuple{2,}, F::NTuple{2,})
    # check quantum number isa halfintrger ?
    all(ishalfinteger, (L..., S..., J..., I..., F...)) || error("transition_matrix_element_m1_total: ")
    L = HalfInt.(L)
    S = HalfInt.(S)
    J = HalfInt.(J)
    I = HalfInt.(I)
    F = HalfInt.(F)
    transition_matrix_element_m1_total(L, S, J, I, F)
end


@doc raw"""

    transition_matrix_element_e1_total(L, S, J, I, F)

compute `c` of ED transition matrix element $\langle \alpha JIF\left\|\mathrm{r}\right\|\alpha' J' I' F' \rangle$.
```math
\langle \alpha JIF\left\|\mathrm{r}\right\|\alpha' J' I' F' \rangle = C * \langle \alpha L\left\|\mathrm{r}\right\|\alpha L'\rangle
```
"""
function transition_matrix_element_e1_total(L::NTuple{2, T}, S::NTuple{2, T}, J::NTuple{2, T}, 
        I::NTuple{2, T}, F::NTuple{2, T}) where {T <: HalfInteger}
    #   <JIF||T^1||J'I'F'> = c2 * <J||T^1||J>
    c2 = uncoup_T1k(J[1], I[1], F[1], J[2], I[2], F[2], 1)
    #   <LSJ||r||L'S'J'> = c3 * <L||r||L>
    #   <L||r||L'> depends on ψ(r), radical wave function.
    c3 = uncoup_T1k(L[1], S[1], J[1], L[2], S[2], J[2], 1)
    return c2 * c3
end
function transition_matrix_element_e1_total(L::NTuple{2,}, S::NTuple{2,}, J::NTuple{2,}, 
    I::NTuple{2,}, F::NTuple{2,})
    # check quantum number isa halfintrger ?
    all(ishalfinteger, (L..., S..., J..., I..., F...)) || error("transition_matrix_element_e1_total: ")
    L = HalfInt.(L)
    S = HalfInt.(S)
    J = HalfInt.(J)
    I = HalfInt.(I)
    F = HalfInt.(F)
    transition_matrix_element_e1_total(L, S, J, I, F)
end