using ZeemanSpectra

using Unitful
using AtomData

using GLMakie
##

atom = Li6

##
B_range = (0:0.01:10)u"Gauss"
E_mat = []
for B in B_range
    Ev = zeeman_struc(atom, 2, B).ES
    if B == B_range[1]
        E_mat = Ev
    else
        E_mat = hcat(E_mat, Ev)
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

lu = 2
ll = 0

if atom == Li6
    if (lu, ll) == (1, 0)
        kx = collect(14903.1:0.00001:14903.5)u"cm^-1"
    elseif (lu, ll) == (2, 0)
        kx = collect(14903.52:0.00001:14903.75)u"cm^-1"
    end
end

T = 1u"K"
P = 0u"Torr"
νp = 2 * 5u"MHz/Torr" * P

df_spec = zeeman_spec(atom, lu, ll, kx, T, νp, BF)

let
    fig = Figure()
    ax = Axis(fig[1, 1])

    kxu = ustrip(kx)
    offset = round(Int, kxu[1])

    line_p = lines!(ax, kxu .- offset, sum(df_spec[df_spec.q .== 0, :c]), label = "q = 0")
    line_s = lines!(ax, kxu .- offset, sum(df_spec[abs.(df_spec.q) .== 1, :c])/2, label = "q = ±1")

    ax.xtickformat = xs -> ["$(x + offset)" for x in xs]
    axislegend()
    fig
end