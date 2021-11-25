using ZeemanSpectra

using Unitful
using AtomData

using GLMakie
##

atom = I127

##
B_range = (0:1:1000)u"Gauss"
E_mat = []
for B in B_range
    E = zeeman_struc(atom, 0, B).ES
    if B == B_range[1]
        E_mat = E
    else
        E_mat = hcat(E_mat, E)
    end
end

let
    fig = Figure()
    ax = Axis(fig[1, 1])

    [lines!(ax, ustrip.(B_range), ustrip.(E)) for E in eachrow(E_mat)]
    
    fig
end

##
BF = 0.0u"Gauss"

kx = collect(7602.2:0.0001:7603.8)u"cm^-1"

T = 200u"K"
P = 0u"Torr"
νp = 2 * 5u"MHz/Torr" * P

df_spec = zeeman_spec(atom, 0, 1, kx, T, νp, BF)

let
    fig = Figure()
    ax = Axis(fig[1, 1])

    kxu = ustrip(kx)
    offset = round(Int, kxu[1])

    cp = sum(df_spec[df_spec.q .== 0, :c])
    cs = sum(df_spec[abs.(df_spec.q) .== 1, :c])/2

    line_p = lines!(ax, kxu .- offset, cp, label = "q = 0")
    line_s = lines!(ax, kxu .- offset, cs, label = "q = ±1")

    ax.xtickformat = xs -> ["$(x + offset)" for x in xs]
    axislegend()
    fig
end

##
L = 1
S = 1/2
J = (3/2, 1/2)
I = 5/2
F = (4, 3)
s1 = (L, S, J[1], I, F[1])
s2 = (L, S, J[2], I, F[2])
k34 = k_I127(F[1], F[2])
A34 = einstein_A_M1(k34, s1, s2)