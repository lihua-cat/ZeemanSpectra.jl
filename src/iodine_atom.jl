struct Iodine <: Atom
    M::Mass                         # relative atomic mass, default in unit `u`
    gI::Float64                     # nuclear Lande g-factor
    I::HalfInt                      # nuclear spin quantum number, either an integer or a half-integer
    ground::FineLevel               # only one ground state 2P3/2
    excited::FineLevel              # only consider first excited states 2P1/2
    function Iodine(M, gI, I, ground, excited)
        if M.val > 0 && I >= 0
            if getproperty.(Ref(ground), [:L, :S, :J]) == [1, 1/2, 3/2]    # 2P3/2 
                if getproperty.(Ref(excited), [:L, :S, :J]) == [1, 1/2, 1/2]    # 2P1/2
                    new(M, gI, I, ground, excited)
                else
                    error("for iodine atoms, first excited state term must be 2P1/2")
                end
            else
                error("for iodine atoms, ground state term must be 2P3/2")
            end
        else
            error("M: $M or I: $I is negative!!")
        end
    end
end

function Iodine(;M, gI, I, A, B, E)
    A = A * u"cm^-1" / 1000
    B = B * u"cm^-1" / 1000
    E = E * u"cm^-1"
    L = [1, 1]
    S = Vector{HalfInt}([1/2, 1/2])
    J = Vector{HalfInt}([3/2, 1/2])
    ground  = FineLevel(L[1], S[1], J[1], A[1], B[1], E[1])
    excited = FineLevel(L[2], S[2], J[2], A[2], B[2], E[2])
    Iodine(M, float(gI), HalfInt(I), ground, excited)
end

I127 = Iodine(M = 126.9u"u", gI = 1.1232, I = 5//2, A = [27.59, 219.73], B = [38.12, 0], E = [0, 7602.9768])