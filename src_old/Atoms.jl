module Atoms

using Unitful, HalfIntegers
import Unitful: Wavenumber, Mass, Time

abstract type Atom end
abstract type EnergyLevel end

# using MyShow
# @add_myshow_type Atom
# @add_myshow_type EnergyLevel

#=
quantum number notation
    L ---- electron orbital quantum number
    S ---- electron spin quantum number
    J ---- L + S, total quantum number of fine structure
    I ---- nuclear spin quantum number
    F ---- J + I, total quantum number of hyperfine sturcture 
=#
struct FineLevel <: EnergyLevel
    L::Int
    S::HalfInt
    J::HalfInt
    A::Wavenumber
    B::Wavenumber
    E::Wavenumber
    function FineLevel(L, S, J, A, B, E)
        if all((L, S, J) .>= 0)
            if J in abs(L-S):abs(L+S)
                new(L, S, J, A, B, E)
            else
                error("L, S, J not match: $L, $S, $J")
            end
        else
            error("L: $L, S: $S, J: $J must be non-negative!!")
        end
    end
end

export Li6
include("alkali_atom.jl")

export I127
include("iodine_atom.jl")

end