module ZeemanSpectra

using Unitful, HalfIntegers, LinearAlgebra
using DataFrames
using WignerSymbols, RationalRoots
import PhysicalConstants.CODATA2018: c_0, h, μ_B, m_e, m_p, μ_0

#   physical constants
const c = ustrip(u"cm/s", c_0)
const me = ustrip(u"u", m_e)
const mp = ustrip(u"u", m_p)
const μB = ustrip(u"Gauss^-1*cm^-1", μ_B/ h/ c_0)
const μN = μB * me / mp

"aa(j::HalfInteger) = Rational(j * (j + 1))"
const aa(j::HalfInteger) = Rational(j * (j + 1))


export ALKALI, IODINE, ATOM
const ALKALI = Set(["Li6", "K39", "K41", "Rb85", "Rb87", "Cs133"])
const IODINE = Set(["I127"])
const ATOM = union(ALKALI, IODINE)


# ------------- import ATOM_DATA ------------------------------------
export ATOM_DATA

include("AtomData.jl")
using .AtomData
# ------------- import LineProfile ---------------------------------
export fwhm_doppler, fwhm_pressure
export profile_voigt, profile_doppler, profile_pressure

using LineProfile
# ------------- angular operator ------------------------------------
export jplus, jminus, angular

include("angular_operator.jl")
# ------------- struct ----------------------------------------------
export zeeman_struc

include("zeeman_struc.jl")
# ------------- WignerEckart ----------------------------------------
export we_Tkq, rme_J, uncoup_T1k, uncoup_T2k

include("wigner_eckart.jl")
#-------------- transition_matrix_element ---------------------------
export transition_matrix_element_m1, transition_matrix_element_e1
export transition_matrix_element_m1_total, transition_matrix_element_e1_total

include("transition_matrix_element.jl")
#-------------- transition_cross_section  ---------------------------
export einstein_A_m1, k_I127, A_I127, σ0_I127, line_I127

include("transition_cross_section.jl")
#-------------- spectra ---------------------------------------------
export zeeman_spec

include("zeeman_spec.jl")
#--------------- END ------------------------------------------------

end
