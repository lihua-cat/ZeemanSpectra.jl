function jplus(j::T, m::T) where {T<:HalfInteger}
    (j > 0 && abs(m) <= j) || error("jplus: j=$j, m=$m")
    RationalRoot(√((j - m) * (j + m + 1)))
end

function jminus(j::T, m::T) where {T<:HalfInteger}
    (j > 0 && abs(m) <= j) || error("jminus: j=$j, m=$m")
    RationalRoot(√((j + m) * (j - m + 1)))
end

function angular(J1::T, I1::T, MJ1::T, MI1::T, 
                 J2::T, I2::T, MJ2::T, MI2::T, 
                 s::Vector{String}) where {T<:HalfInteger}
    (J1 == J2 && I1 == I2) || error("angular: J1=$(J1), J2=$(J2), I1=$(I1), I1=$(I1)")
    e = 1
    for op in reverse(s)
        if op == "Jz"
            e *= MJ2
        elseif op == "Iz"
            e *= MI2
        elseif op == "Jp"
            e *= jplus(J2, MJ2)
            MJ2 += 1
        elseif op == "Ip"
            e *= jplus(I2, MI2)
            MI2 += 1
        elseif op == "Jm"
            e *= jminus(J2, MJ2)
            MJ2 -= 1
        elseif op == "Im"
            e *= jminus(I2, MI2)
            MI2 -= 1
        else
            error("no such operator!")
        end
        if abs(MJ2) > J2 || abs(MI2) > I2
            return 0
        end
    end 
    return (MJ1, MI1) == (MJ2, MI2) ? e : 0
end