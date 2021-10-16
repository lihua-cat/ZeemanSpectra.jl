module AtomData

using Unitful, DataFrames, HalfIntegers

export ATOM_DATA

#   isotope names
Name = ["Li6", "K39", "K41", "Rb85", "Rb87", "Cs133", "I127"]

#   relative atomic mass(1 u=1.66053906660(50)e−24 g)
M = [6.0151214, 38.9637, 40.9618, 84.9118, 86.9092, 132.9055, 126.9]u"u"

#   nuclear spin quantum number
I = HalfInt.([1, 3 // 2, 3 // 2, 5 // 2, 3 // 2, 7 // 2, 5 // 2])

#   nuclear Lande g-factor
gI = [0, 0.261, 0.143, 0.539, 1.827, 0.732, 1.1232]

#   D1 line wavelength, cm-1
D1 = [14903.298, 12985.1851928, 12985.193050, 12578.9483900, 12578.95098147, 11178.26815870, 7602.9768]u"cm^-1"

#   D2 line wavelength, cm-1
D2 = [14903.633, 13042.8954964, 13042.903375, 12816.54678496, 12816.54938993, 11732.3071041, missing]u"cm^-1"

#   hyperfine split constant A in ground state, cm-1
Ag = [5.1, 7.700660702, 4.236597573, 33.75368436, 113.9901925, 76.6583661, 27.59]u"cm^-1" / 1000

#   hyperfine split constant B in ground state(0 if J=1/2)
Bg = [0, 0, 0, 0, 0, 0, 38.12]u"cm^-1" / 1000

#   hyperfine split constant A in excited state 1
Ae1 = [0.57993, 0.962332415, 0.507017425, 4.026785757, 13.54937355, 9.736735939, 219.73]u"cm^-1" / 1000

#   hyperfine split constant B in excited state 1
Be1 = [0, 0, 0, 0, 0, 0, 0,]u"cm^-1" / 1000

#   hyperfine split constant A in excited state 2
Ae2 = [-0.051702, 0.202139842, 0.113411792, 0.834243802, 2.830291348, 1.679161655, missing]u"cm^-1" / 1000

#   hyperfine split constant B in excited state 2
Be2 = [-0.0033, 0.094398639, 0.111410408, 0.863263878, 0.417622247, -0.012675436, missing]u"cm^-1" / 1000


# dataframe construct
const ATOM_DATA = DataFrame(Name = Name,
                          M    = M,
                          I    = I,
                          gI   = gI,
                          D1   = D1,
                          D2   = D2,
                          Ag   = Ag,
                          Bg   = Bg,
                          Ae1  = Ae1,
                          Be1  = Be1,
                          Ae2  = Ae2,
                          Be2  = Be2)

end # of module