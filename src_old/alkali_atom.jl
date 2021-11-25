struct Alkali <: Atom
    M::Mass                         # relative atomic mass, default in unit `u`
    gI::Float64                     # nuclear Lande g-factor
    I::HalfInt                      # nuclear spin quantum number, either an integer or a half-integer
    ground::FineLevel               # only one ground state 2S1/2
    excited::Vector{FineLevel}      # several excited states
    function Alkali(M, gI, I, ground, excited)
        if M.val > 0 && I >= 0
            if getproperty.(Ref(ground), [:L, :S, :J]) == [0, 1/2, 1/2]  # 2S1/2
                if all(getproperty.(excited, :S) .== 1/2)
                    new(M, gI, I, ground, excited)
                else
                    error("for alkali atoms, S: $S must be 1/2")
                end
            else
                error("for alkali atoms, ground state term must be 2S1/2")
            end
        else
            error("M: $M or I: $I is negative!!")
        end
    end
end

function Alkali_D1_D2(;M, gI, I, A, B, E)
    A = A * u"cm^-1" / 1000
    B = B * u"cm^-1" / 1000
    E = E * u"cm^-1"
    L = [0, 1, 1]
    S = Vector{HalfInt}([1/2, 1/2, 1/2])
    J = Vector{HalfInt}([1/2, 1/2, 3/2])
    ground  = FineLevel(L[1], S[1], J[1], A[1], B[1], E[1])
    excited = [FineLevel(L[2], S[2], J[2], A[2], B[2], E[2]), FineLevel(L[3], S[3], J[3], A[3], B[3], E[3])]
    Alkali(M, float(gI), HalfInt(I), ground, excited)
end

Li6 = Alkali_D1_D2(M = 6.0151214u"u", gI = 0, I = 1, A = [5.1, 0.57993, -0.051702], B = [0, 0, -0.0033], E = [0, 14903.298, 14903.633])
