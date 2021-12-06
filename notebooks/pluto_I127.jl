### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 212e59d1-140e-41f9-bc67-493d2d61d008
begin
	import Pkg
	Pkg.activate("..")
	Base.Text(sprint(io->Pkg.status(io=io)))
end

# ╔═╡ 37294fc5-c716-4e79-81a3-7a4b08cbe599
push!(LOAD_PATH, "C:\\Users\\lihua\\.julia\\dev")

# ╔═╡ 2e60614f-4af7-43bf-b926-1b5425aec7bd
using Revise

# ╔═╡ f1c7f04c-65e5-4cae-b3fe-42e3fefd058d
using ZeemanSpectra, AtomBase

# ╔═╡ e51ad0a1-60b6-4ac0-8a36-853cdf1b036a
using Unitful

# ╔═╡ 811f6a70-a9d6-41bb-bf39-d2a13c26b4cc
using DataFrames

# ╔═╡ 2f3c19bd-7eb7-465a-87ca-2a405c3e7a5d
using Printf

# ╔═╡ b6957132-8575-4a56-b131-9ac3b2d72498
using CairoMakie

# ╔═╡ b1464c4d-47a4-47c6-9676-ba4627260845
html"""<style>
main {
    max-width: 900px;
    align-self: flex-start;
    margin-left: 50px;
}
"""

# ╔═╡ 2a810db4-e771-432d-8880-45667807d507
import PhysicalConstants.CODATA2018: h, μ_B, μ_0

# ╔═╡ faa31604-43ba-4992-bd15-9dc09990ccd3
import PlutoUI

# ╔═╡ 08770a4f-9b10-401f-9bb0-0bd856b67290
CairoMakie.activate!(type = "svg")

# ╔═╡ 57387ba5-2926-4a61-b2fe-26ba61453ddf
md"# Zeeman spectra of $^{127}\mathrm{I}$"

# ╔═╡ f8962e80-eb90-4cf7-9863-458c471a0c39
PlutoUI.TableOfContents(depth = 4)

# ╔═╡ 8d1d5ad5-5f2b-4ffc-8cc5-959aa79ba342
md"## 1. Atom Data"

# ╔═╡ 5f8a9c94-b0b1-4bb3-ae81-a664863dbe53
atom = I127

# ╔═╡ a2bafc87-52ef-4667-b907-9c9d920dd258
md"### Fine Structure"

# ╔═╡ 61054b00-6b53-4fab-abd4-f454f80f59e1
function term(s)
	(;L, S, J) = s
	LL = ("S", "P", "D", "F")
	St = string(2S+1)
	Lt = LL[L+1]
	Jt = string(J)
	t = St*Lt*Jt
end

# ╔═╡ 773bab65-0dcd-4400-9b7f-8dec6effad8e
let
	s0 = atom.ground
	s1 = atom.excited[1]
	A1 = repr(s0.hfc.A, context=:compact=>true)
	B1 = repr(s0.hfc.B, context=:compact=>true)
	A2 = repr(s1.hfc.A, context=:compact=>true)
	B2 = repr(s1.hfc.B, context=:compact=>true)
	md"""
	| Term | Energy | $A_{hf}$ | $B_{hf}$ |
	| ---- | ------ | -------- | -------- |
	|$(term(s0.state))|$(s0.E)|$(A1)|$(B1)|
	|$(term(s1.state))|$(s1.E)|$(A2)|$(B2)|

	All data from Table 2 in [1].
	"""
end

# ╔═╡ e0b645bd-7909-422e-8da1-b6c340926a8f
md"### Hyperfine structure"

# ╔═╡ 3cfa4727-37df-41c8-95fb-bcbd62fa047e
md"ground state: $(term(atom.ground.state))"

# ╔═╡ 537cf404-25f3-49b5-a9b6-09f68384e275
let
	BF = 0.0u"Gauss"
	df0 = zeeman_struc(atom, 0)

	df = DataFrame(F=[], ES = [])
	for F in Set(df0.F)
		ES = filter(row -> row.F == F, df0).ES
		@assert all(ES .≈ ES[1]) "error"
		push!(df, (F, ES[1]))
	end
	sort!(df, :ES, rev=true)
end

# ╔═╡ abbd971e-f498-4f1f-bc88-f6064c3eee64
md"first excited state: $(term(atom.excited[1].state))"

# ╔═╡ 430de0d6-442d-49c1-a2d7-6b22b8b99afc
let
	BF = 0.0u"Gauss"
	df0 = zeeman_struc(atom, 1)

	df = DataFrame(F=[], ES = [])
	for F in Set(df0.F)
		ES = filter(row -> row.F == F, df0).ES
		@assert all(ES .≈ ES[1]) "error"
		push!(df, (F, ES[1]))
	end
	sort!(df, :ES, rev=true)
end

# ╔═╡ 8ade433b-9039-4ef6-9a40-5f219e03b93a
md"see fig.1 in [2]"

# ╔═╡ 2cbf8cb2-95c8-4ed3-b82f-f59e40693ab8
md"## 2. Transition lines"

# ╔═╡ 20a4e829-2847-437e-80db-380c4a728bea
md"### zero field"

# ╔═╡ 48777912-33ef-474d-bcdc-13caac6bb977
md"#### wavenumber, relative intensity and Einstein A coefficient"

# ╔═╡ e348d5ba-29a1-4f53-9608-d47d0afad7df
let
	B = 0.0u"Gauss"
	s1 = atom.ground
	s2 = atom.excited[1]
	L1, L2 = s1.state.L, s2.state.L
	S1, S2 = s1.state.S, s2.state.S
	J1, J2 = s1.state.J, s2.state.J
	I = atom.I
	c2 = reducedME_M1(L1, S1, J1, L2, S2, J2)
	df = DataFrame(F2=[], F1=[], k = [], c = [], A = [])
	for F1 in 4:-1:1, F2 in 3:-1:2
		E1 = filter(row -> row.F == F1, zeeman_struc(atom, 0, B)).E
		E2 = filter(row -> row.F == F2, zeeman_struc(atom, 1, B)).E
		@assert all(E1 .≈ E1[1]) "error 1"
		@assert all(E2 .≈ E2[1]) "error 2"
		k = E2[1] - E1[1]
		c1 = uncoup_T1(J1, I, F1, J2, I, F2, 1)
		c = (c1 * c2)^2
		c1 ≈ 0 && continue
		A = aᵢⱼ(k, c) / (2F2 + 1)
		push!(df, (F1, F2, k, c1^2, A))
	end
	df.c = df.c / sum(df.c)
	sort!(df, [:c, :A], rev=true)
	df.k = [@sprintf("%.4f", k.val) for k in df.k]
	df
end

# ╔═╡ 5b049fac-21f3-4053-8d09-da7ae5a8fd22
begin
	A = aᵢⱼ(atom.excited[1].E, reducedME_M1(1, 1/2, 3/2, 1, 1/2, 1/2)^2) / (2*1/2+1)
	md"""
	``{}^2P_{1/2} \rightarrow {}^2P_{3/2}``: $(A)

	A in [3]: 7.87 s^-1
	
	k and c: see Table 4 in [1]
	"""
end

# ╔═╡ df5f0239-6f1a-4e6b-9236-9738dbc20fdd
uncoup_T1(3/2, 5/2, 4, 1/2, 5/2, 3, 1)^2 * 7.93979399211916 * 2 / 7

# ╔═╡ 15417b8a-fd8c-4174-9bea-7aadd422660f
uncoup_T1(3/2, 5/2, 2, 1/2, 5/2, 2, 1)^2 * 7.939805013011782 * 2 / 5

# ╔═╡ 70d1790d-7bb9-456d-82a3-cadae892803c
md"#### transition cross section"

# ╔═╡ f5956f80-ff23-434d-9d81-64aee5a6032b
γ = let
	γ_He = 4.1u"MHz/Torr"
	γ_O₂ = 5.7u"MHz/Torr"
	(γ_He * 0.8 + γ_O₂ * 0.2)
end

# ╔═╡ 4155d55a-ea9c-4343-ae56-69705b94984b
let
	BF = 0.0u"Gauss"
	kx = collect(7602.2:0.0001:7603.8)u"cm^-1"
	
	T = 200u"K"
	P = 10u"Torr"
	νp = 2 * γ * P
	
	df_spec = zeeman_spec(atom, 0, 1, kx, T, νp, "M1", BF)
	df_new = DataFrame(F1=[], F2=[], k=[], c=[], A=[], σ0=[], σ=[])

	for F1 in Set(df_spec.F1), F2 in Set(df_spec.F2)
		df = filter(row->row.F1==F1&&row.F2==F2, df_spec)
		c = sum(df.c0)
		c < eps() && continue
		@assert all(df.k0 .≈ df.k0[1]) "error"
		A = sum(df.a)/(2F2+1)
		σ0 = sum(df.σ0)/(2F2+1)
		σ = sum(df.σ)/(2F2+1)
		push!(df_new, (F1, F2, df.k0[1], c, A, σ0, σ))
	end
	sort!(df_new, [:c, :A], rev=true)
	df_new
end

# ╔═╡ cd826b82-1249-4e88-ac2b-7a25bfe24f61
let
	σ = 1.3e-16/√(200)
	md"""
	in [4], ``F'=3 \rightarrow F=4, \sigma = 1.3 \times 10^{16}/\sqrt{T}``(has been multiplied by 7/12)
	
	T = 200K, $\sigma/(7/12) =$ $(σ/(7/12))
	"""
end

# ╔═╡ 5112acd1-614a-48d3-b826-9a6a56b3ade9
σ0_I127(4, 3, T = 200u"K", P = 0u"Torr")

# ╔═╡ a524ef95-3382-4b14-bc5a-6a4d79c39c67
line_I127(4, 3)

# ╔═╡ 769af577-cdb9-4780-9446-1e45c6c45936
md"### nonzero field"

# ╔═╡ 579d0e52-fe79-4678-a21e-326b6773d341
md"#### Energy splits"

# ╔═╡ 55440c43-8999-484d-bbf5-3d234841289b
let
	B_range = (0:50:4000)u"Gauss"
	E_mat = []
	for B in B_range
	    E = zeeman_struc(atom, 0, B).ES
	    if B == B_range[1]
	        E_mat = E
	    else
	        E_mat = hcat(E_mat, E)
	    end
	end
	F_list = zeeman_struc(atom, 0, B_range[1]).F
	
	fig = Figure()
    ax = Axis(fig[1, 1])

    color = [:blue, :orange, :green, :red]

    for n in 1:size(E_mat, 1)
        lines!(ax, ustrip.(B_range), ustrip.(E_mat[n, :]), label = "$(F_list[n])", color = color[Int(F_list[n])])
    end

	elem = [LineElement(color = color[i]) for i in 4:-1:1]
	label = ["$(unique(F_list)[i])" for i in 1:4]
	
    Legend(fig[1, 2], elem, label, "F")
    fig
end

# ╔═╡ a2f787ed-ed51-441d-a58c-b342cdfb7a29
md"#### relative transition cross-section"

# ╔═╡ 22251140-ae5d-4736-a926-6729f6d9ea29
T = 200u"K"

# ╔═╡ 8cba09c5-4a1c-4bc5-aada-2df6da6a52a5
P = 10u"Torr"

# ╔═╡ 4fb83bab-da18-4934-a7ab-000979952ad9
let
	BF = 0u"Gauss"

	kx = collect(7602.2:0.0001:7603.8)u"cm^-1"
	
	νp = 2 * γ * P
	
	df0 = zeeman_spec(atom, 0, 1, kx, T, νp, "M1", zero(BF))
	
	σ34 = σ0_I127(4, 3, T=T, P=P, γ=γ) * 7 / 3
	
	df_spec = zeeman_spec(atom, 0, 1, kx, T, νp, "M1", BF)

    fig = Figure(resolution = (1200, 600), fontsize = 24, font = "sans")
    ax1 = Axis(fig[1, 1], title = "P", palette = (color = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#e377c2", "#bcbd22", "#17becf"],))
	ax2 = Axis(fig[1, 2], title = "S", palette = (color = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#e377c2", "#bcbd22", "#17becf"],))
	hideydecorations!(ax2; label = true, ticklabels = true, ticks = true, grid = false, minorgrid = false, minorticks = false)
    linkyaxes!(ax1, ax2)

    df_p = filter(row->row.q == 0, df_spec)
	df_s = filter(row->abs(row.q) == 1, df_spec)

    kxu = ustrip(kx)
    offset = round(Int, kxu[1])

    for F1 in Set(df_p.F1), F2 in Set(df_p.F2)
		p = sum(filter(row->row.F1==F1&&row.F2==F2, df_p).σ)/σ34
        s = sum(filter(row->row.F1==F1&&row.F2==F2, df_s).σ)/σ34/2
        lines!(ax1, kxu .- offset, p, label = L"%$F2 \rightarrow %$F1")
		lines!(ax2, kxu .- offset, s, label = L"%$F2 \rightarrow %$F1")
    end
	
    ax1.xtickformat = xs -> ["$(x + offset)" for x in xs]
	ax2.xtickformat = xs -> ["$(x + offset)" for x in xs]
	ax1.yticks = 0:0.2:1
    axislegend(ax2)
	ylims!(-0.02, 1.02)
	
    fig
end

# ╔═╡ a9ae9150-9ad1-44fd-9773-360710f1c4f5
let
	B_range = collect(0:10:600)u"Gauss"
	F1_list = 4:-1:1
	F2_list = 3:-1:2
	kx = collect(7602.2:0.0001:7603.8)u"cm^-1"
	σ34_0 = σ0_I127(4, 3, T = T, P = P, γ = γ)
	σ_ma = Matrix{Float64}(undef, length(B_range), length(F1_list)*length(F2_list))
	for i in 1:length(B_range)
		df = zeeman_spec(atom, 0, 1, kx, T, 2P*γ, "M1", B_range[i])
	    j = 1
	    for F1 in F1_list
	        for F2 in F2_list
	            σ = sum(filter(row -> row.F1==F1&&row.F2==F2&&abs(row.q)==1, df).σ)/2*3
	            σm = maximum(σ)
	            σ_ma[i, j] = σm/(σ34_0 * 7)
	            j += 1
	        end
	    end
	end

    fig = Figure(resolution = (1200, 600), fontsize = 24, font = "sans")

    ax = Axis(fig[1,1], xlabel = "Magnetic B field(Gauss)", ylabel = "Relative SSG", palette = (color = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#e377c2", "#bcbd22", "#17becf"],))

    i = 1
    Bx = ustrip.(u"Gauss", B_range)
	for F1 in F1_list, F2 in F2_list
		σm = ustrip.(NoUnits, σ_ma[:, i])
		lines!(ax, Bx, σm, linewidth = 2, label = L"%$F2\rightarrow%$F1")
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
    # save("./plot/relatvie_ssg_I127.svg", fig, px_per_unit = 0.5)
    fig
end

# ╔═╡ f4f8fc03-2ce6-4bd1-8bce-f7f04fe12854
σm_I127(3, 3, 0u"Gauss", "P"; T = T, P = P, γ = γ) / 7

# ╔═╡ bc9aa2fd-e9af-4f85-a7d7-587810ef638e
σm_I127(3, 3, 0u"Gauss", "S"; T = T, P = P, γ = γ) / 7

# ╔═╡ f33dfb79-0b28-4d14-95b5-4c49c075114c
md"""
## Reference
1. Engleman, R., Keller, R. A. & Palmer, B. A. Hyperfine structure and isotope shift of the 1.3-μm transition of ^129l. Appl Optics 19, 2767 (1980).
2. Kelly, M. A., McIver, J. K., Shea, R. F. & Hager, G. D. Frequency tuning of a CW atomic iodine laser via the Zeeman effect. Ieee J Quantum Elect 27, 263–273 (1991).
3. Lilenfeld, H. V., Richardson, R. J. & Hovis, F. E. The electron spin resonance spectrum of the 2 P 1/2 state of atomic iodine. J Chem Phys 74, 2129–2132 (1981).
4. Hager, G. D. et al. A simplified analytic model for gain saturation and power extraction in the flowing chemical oxygen-iodine laser. Ieee J Quantum Elect 32, 1525–1536 (1996).
  
"""

# ╔═╡ Cell order:
# ╠═b1464c4d-47a4-47c6-9676-ba4627260845
# ╠═2e60614f-4af7-43bf-b926-1b5425aec7bd
# ╠═37294fc5-c716-4e79-81a3-7a4b08cbe599
# ╠═212e59d1-140e-41f9-bc67-493d2d61d008
# ╠═f1c7f04c-65e5-4cae-b3fe-42e3fefd058d
# ╠═e51ad0a1-60b6-4ac0-8a36-853cdf1b036a
# ╠═2a810db4-e771-432d-8880-45667807d507
# ╠═811f6a70-a9d6-41bb-bf39-d2a13c26b4cc
# ╠═2f3c19bd-7eb7-465a-87ca-2a405c3e7a5d
# ╠═faa31604-43ba-4992-bd15-9dc09990ccd3
# ╠═b6957132-8575-4a56-b131-9ac3b2d72498
# ╠═08770a4f-9b10-401f-9bb0-0bd856b67290
# ╟─57387ba5-2926-4a61-b2fe-26ba61453ddf
# ╟─f8962e80-eb90-4cf7-9863-458c471a0c39
# ╟─8d1d5ad5-5f2b-4ffc-8cc5-959aa79ba342
# ╠═5f8a9c94-b0b1-4bb3-ae81-a664863dbe53
# ╟─a2bafc87-52ef-4667-b907-9c9d920dd258
# ╟─61054b00-6b53-4fab-abd4-f454f80f59e1
# ╟─773bab65-0dcd-4400-9b7f-8dec6effad8e
# ╟─e0b645bd-7909-422e-8da1-b6c340926a8f
# ╟─3cfa4727-37df-41c8-95fb-bcbd62fa047e
# ╟─537cf404-25f3-49b5-a9b6-09f68384e275
# ╟─abbd971e-f498-4f1f-bc88-f6064c3eee64
# ╟─430de0d6-442d-49c1-a2d7-6b22b8b99afc
# ╟─8ade433b-9039-4ef6-9a40-5f219e03b93a
# ╟─2cbf8cb2-95c8-4ed3-b82f-f59e40693ab8
# ╟─20a4e829-2847-437e-80db-380c4a728bea
# ╟─48777912-33ef-474d-bcdc-13caac6bb977
# ╟─e348d5ba-29a1-4f53-9608-d47d0afad7df
# ╟─5b049fac-21f3-4053-8d09-da7ae5a8fd22
# ╠═df5f0239-6f1a-4e6b-9236-9738dbc20fdd
# ╠═15417b8a-fd8c-4174-9bea-7aadd422660f
# ╟─70d1790d-7bb9-456d-82a3-cadae892803c
# ╟─f5956f80-ff23-434d-9d81-64aee5a6032b
# ╟─4155d55a-ea9c-4343-ae56-69705b94984b
# ╟─cd826b82-1249-4e88-ac2b-7a25bfe24f61
# ╠═5112acd1-614a-48d3-b826-9a6a56b3ade9
# ╠═a524ef95-3382-4b14-bc5a-6a4d79c39c67
# ╟─769af577-cdb9-4780-9446-1e45c6c45936
# ╟─579d0e52-fe79-4678-a21e-326b6773d341
# ╟─55440c43-8999-484d-bbf5-3d234841289b
# ╟─a2f787ed-ed51-441d-a58c-b342cdfb7a29
# ╟─22251140-ae5d-4736-a926-6729f6d9ea29
# ╟─8cba09c5-4a1c-4bc5-aada-2df6da6a52a5
# ╟─4fb83bab-da18-4934-a7ab-000979952ad9
# ╠═a9ae9150-9ad1-44fd-9773-360710f1c4f5
# ╠═f4f8fc03-2ce6-4bd1-8bce-f7f04fe12854
# ╠═bc9aa2fd-e9af-4f85-a7d7-587810ef638e
# ╟─f33dfb79-0b28-4d14-95b5-4c49c075114c
