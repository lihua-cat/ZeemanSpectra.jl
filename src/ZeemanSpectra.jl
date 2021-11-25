module ZeemanSpectra

import PhysicalConstants.CODATA2018: c_0 as 𝑐, h as ℎ, μ_B, m_e as 𝑚e, m_p as 𝑚p, μ_0 as 𝜇0
using Unitful
import Unitful: Wavenumber, Frequency, Area, BField, Quantity
using HalfIntegers
using DataFrames

import Base: nameof
nameof(::Type{Wavenumber}) = :Wavenumber
nameof(::Type{Frequency}) = :Frequency
nameof(::Type{Area}) = :Area

using AtomBase
using AtomData
using LineProfile

const 𝜇B = μ_B / (ℎ * 𝑐) |> u"cm^-1/Gauss"
const 𝜇N = 𝜇B * 𝑚e / 𝑚p
const gS = 1.00115965218085 * 2 #https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.030801

include("hamiltonian.jl")

export zeeman_struc
include("zeeman structure.jl")

export zeeman_spec
include("zeeman spectra.jl")

export einstein_A_M1, k_I127, A_I127, σ0_I127, line_I127
include("transition cross section.jl")

end