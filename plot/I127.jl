##
using ZeemanSpectra
using Unitful
using CairoMakie

atom = I127

T = 200u"K"
P = 10u"Torr"

γ_He = 4.1u"MHz/Torr"
γ_O₂ = 5.7u"MHz/Torr"
γ = (γ_He * 0.8 + γ_O₂ * 0.2)

σ34_0 = σ0_I127(4, 3, T = T, P = P, γ = γ)

spec_I127(kx, T, P, BF) = zeeman_spec(atom, 0, 1, kx, T, 2P*γ, "M1", BF)
df_F1_F2(F1, F2, df) = filter(row -> row.F1 == F1&&row.F2==F2, df)

##
F1, F2 = 4, 3
k0 = k_I127(F1, F2)
kx_0 = collect(k0-0.05u"cm^-1":0.0001u"cm^-1":k0+0.05u"cm^-1")
BF = 400u"Gauss"
df_spec = spec_I127(kx_0, T, P, BF)
df = df_F1_F2(F1, F2, df_spec)

##
let 
    cycle = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    offset = round(ustrip.(u"cm^-1", k0), digits = 2)
    x = ustrip.(u"cm^-1", kx_0) .- offset
    fig = Figure(resolution = (1200, 600), fontsize = 24, font = "sans")

    ax1 = Axis(fig[1,1], xlabel = L"\Delta k(cm^{-1})", ylabel = "Relative transition cross-section", xticks = LinearTicks(6))
    ax2 = Axis(fig[1,2], xlabel = L"\Delta k(cm^{-1})", ylabel = "", xticks = LinearTicks(6))
    hideydecorations!(ax2; label = true, ticklabels = true, ticks = true, grid = false,
        minorgrid = false, minorticks = false)
    linkyaxes!(ax1, ax2)

    σp = zeros(length(kx_0))
    σs = zeros(length(kx_0))
    for row in eachrow(df)
        k = ustrip.(u"cm^-1", row.k0) - offset
        if row.q == 0
            σ = ustrip.(NoUnits, row.σ*3/σ34_0/(2F2+1))
            σm = maximum(σ)
            lc = cycle[1]
            σp .+= σ
            lines!(ax1, x, σ, lw = 1, color = lc)
            # text!(ax1, L"%$(row.MF2)\rightarrow%$(row.MF1)", position = (k, σm), textsize = 16, align = (:center, :baseline))
        else
            σ = ustrip.(NoUnits, row.σ/2*3/σ34_0/(2F2+1))
            σm = maximum(σ)
            σs .+= σ
            if row.q == 1
                lc = cycle[2]
            elseif row.q == -1
                lc = cycle[3]
            end
            lines!(ax2, x, σ, lw = 1, color = lc)
            # text!(ax2, L"%$(row.MF2)\rightarrow%$(row.MF1)", position = (k, σm), textsize = 16, align = (:center, :baseline))
        end
    end
    lines!(ax1, x, σp, linewidth = 1, color = cycle[4], linestyle = :dash)
    lines!(ax2, x, σs, linewidth = 1, color = cycle[5], linestyle = :dash)
    # xlims!(ax1, x[1], x[end])
    ylims!(ax1, 0, nothing)

    Label(fig[1, 1], "(a)", textsize = 30, tellheight = false, tellwidth = false, halign = :left, valign = :top, padding = (10, 10, 10, 10))
    Label(fig[1, 2], "(b)", textsize = 30, tellheight = false, tellwidth = false, halign = :left, valign = :top, padding = (10, 10, 10, 10))

    elem_1 = LineElement(color = cycle[1])
    elem_2 = LineElement(color = cycle[2])
    elem_3 = LineElement(color = cycle[3])
    elem_4 = LineElement(color = cycle[4], linestyle = :dash)
    elem_5 = LineElement(color = cycle[5], linestyle = :dash)
    Legend(fig[1, 1], [elem_1, elem_4], [L"\Delta M_F = 0", L"P"], tellwidth = false, halign = :right, valign = :top, labelhalign = :center, framevisible = true, labelsize = 20)
    Legend(fig[1, 2], [elem_2, elem_3, elem_5], [L"\Delta M_F = +1", L"\Delta M_F = -1", L"S"], tellwidth = false, halign = :right, valign = :top, labelhalign = :center, framevisible = true, labelsize = 20)
    # save("./plot/relative_cross_section_I127.svg", fig, px_per_unit = 0.5)
    # save("./plot/relative_cross_section_I127.png", fig, px_per_unit = 0.5)
    fig
end

##
B_range = collect(0:1:600)u"Gauss"
F1_list = 4:-1:1
F2_list = 3:-1:2
kx = collect(7602.2:0.0001:7603.8)u"cm^-1"
σ_ma = Matrix{Float64}(undef, length(B_range), length(F1_list)*length(F2_list))
for i in 1:length(B_range)
    df = spec_I127(kx, T, P, B_range[i])
    j = 1
    for F1 in F1_list
        for F2 in F2_list
            σ = sum(filter(row->abs(row.q) == 1, df_F1_F2(F1, F2, df)).σ)
            σm = maximum(σ)
            σ_ma[i, j] = σm/σ34_0/7/2*3
            j += 1
        end
    end
end

##
let
    cycle = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#bcbd22", "#17becf"]
    fig = Figure(resolution = (1200, 600), fontsize = 24, font = "sans")

    ax = Axis(fig[1,1], xlabel = "Magnetic B field(Gauss)", ylabel = "Relative SSG", xticks = LinearTicks(6))

    i = 1
    Bx = ustrip.(u"Gauss", B_range)
	for F1 in F1_list, F2 in F2_list
		σm = ustrip.(NoUnits, σ_ma[:, i])
		lines!(ax, Bx, σm, linewidth = 2, color = cycle[i], label = L"%$F2\rightarrow%$F1")
		i += 1
	end
    # gth_to_g0 = 0.34
    # hlines!(ax3, [gth_to_g0], color = (:blue, 0.5), linestyle = [10, 16, 24])
    # pos = (Bx[is], σ_ma[is, 1])
    # scatter!(ax, pos, markersize = 8, color = :red)
    # text!(ax, "($(Bx[is]), $(round(Float64(σ_ma[is, 1]), digits = 4)))", position = pos.*1.05, textsize = 18, color = :red)
    fig[1, 1] = Legend(fig, ax, L"F' \rightarrow F", tellwidth = false, halign = :right, valign = :top)
    xlims!(-10, 610)
    # ylims!(-0.02, 1)
    ax.yticks = 0:0.2:1.0
    ax.xticks = 0:50:600
    save("./plot/relatvie_ssg_I127.svg", fig, px_per_unit = 0.5)
    fig
end