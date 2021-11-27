using ZeemanSpectra

using Unitful
using AtomData

using GLMakie
##

atom = Li6

##
B_range = (0:0.01:10)u"Gauss"
level = 2
E_mat = []
for B in B_range
    Ev = zeeman_struc(atom, level, B).ES
    if B == B_range[1]
        E_mat = Ev
    else
        E_mat = hcat(E_mat, Ev)
    end
end

F_list = zeeman_struc(atom, level, B_range[1]).F

let
    fig = Figure()
    ax = Axis(fig[1, 1])

    color = [:blue, :orange, :green, :red]

    for n in 1:size(E_mat, 1)
        lines!(ax, ustrip.(B_range), ustrip.(E_mat[n, :]), label = "$(F_list[n])", color = color[Int(F_list[n]+1/2)])
    end
    axislegend(ax)
    fig
end

##
BF = 0.0u"Gauss"

lu = 1
ll = 0

if atom == Li6
    if (lu, ll) == (1, 0)
        kx = collect(14903.1:0.00001:14903.5)u"cm^-1"
    elseif (lu, ll) == (2, 0)
        kx = collect(14903.52:0.00001:14903.75)u"cm^-1"
    end
end

T = 50u"K"
P = 0u"Torr"
νp = 2 * 5u"MHz/Torr" * P

df_spec = zeeman_spec(atom, ll, lu, kx, T, νp, "E1", BF)

let
    fig = Figure()
    ax = Axis(fig[1, 1])

    kxu = ustrip(kx)
    offset = round(Int, kxu[1])

    for F1 in Set(df_spec.F1), F2 in Set(df_spec.F2)
        p = sum(filter(row->row.F1==F1&&row.F2==F2, df_spec).c)
        lines!(ax, kxu .- offset, p, label = "$F2 -> $F1")
    end
    # line_p = lines!(ax, kxu .- offset, sum(df_spec[df_spec.q .== 0, :c]), label = "q = 0")
    # line_s = lines!(ax, kxu .- offset, sum(df_spec[abs.(df_spec.q) .== 1, :c])/2, label = "q = ±1")

    ax.xtickformat = xs -> ["$(x + offset)" for x in xs]
    axislegend()
    fig
end