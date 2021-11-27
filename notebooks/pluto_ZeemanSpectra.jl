### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 49bb3630-4e0d-11ec-39f2-c3d1023de02a
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add(url = "https://github.com/lihua-cat/UsefulFunctions.jl")
	Pkg.add(url = "https://github.com/lihua-cat/AtomData.jl")
	Pkg.add(url = "https://github.com/lihua-cat/AtomBase.jl")
	Pkg.add(url = "https://github.com/lihua-cat/LineProfile.jl")
	Pkg.add(url = "https://github.com/lihua-cat/ZeemanSpectra.jl")
	Pkg.add("DataFrames")
	Pkg.add("Unitful")
	Pkg.add("PhysicalConstants")
end

# ╔═╡ f1c7f04c-65e5-4cae-b3fe-42e3fefd058d
using ZeemanSpectra

# ╔═╡ d1651f70-47c1-4ce7-9f6c-bac9693d12aa
using AtomData, AtomBase

# ╔═╡ e51ad0a1-60b6-4ac0-8a36-853cdf1b036a
using Unitful

# ╔═╡ 811f6a70-a9d6-41bb-bf39-d2a13c26b4cc
using DataFrames

# ╔═╡ 2a810db4-e771-432d-8880-45667807d507
import PhysicalConstants.CODATA2018: h, μ_B, μ_0

# ╔═╡ 5f8a9c94-b0b1-4bb3-ae81-a664863dbe53
atom = I127

# ╔═╡ 1e6a4ab2-837c-4aee-a432-0fec601b126f
atom.ground.hfc

# ╔═╡ b1d4e2b4-d9ad-4ec7-91f7-39080b8bdd4b
atom.excited[1].hfc

# ╔═╡ fe4b4a43-1c6f-4b6e-b21e-0938b41fa657
BF = 0.0u"Gauss"

# ╔═╡ 4155d55a-ea9c-4343-ae56-69705b94984b
struc = zeeman_struc(atom, 0, BF)

# ╔═╡ 1bc2e2c3-f7ef-44e8-832c-b239866f3c11
filter(row -> row.F == 2, struc)

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
		A = 64π^4 / 3h * k^3 * c / (2F2 + 1) |> u"s^-1"
		push!(df, (F2, F1, k, c1^2, A))
	end
	df.c = df.c / sum(df.c)
	df
end

# ╔═╡ ba00a0fa-a7de-493e-a86d-67d132281693
let
	BF = 0.0u"Gauss"

	kx = collect(7602.2:0.0001:7603.8)u"cm^-1"
	
	T = 200u"K"
	P = 10u"Torr"
	νp = 2 * 5u"MHz/Torr" * P
	
	df_spec = zeeman_spec(atom, 1, 0, kx, T, νp, "M1", BF)
	df34 = filter(row->row.F1 == 3 && row.F2 == 4, df_spec)
	df1 = filter(row->row.q == 0, df34)
	df2 = filter(row->row.q == -1, df34)
	df3 = filter(row->row.q == +1, df34)
	[sum(dfi.a) / (2*3+1) for dfi in (df1, df2, df3)]
end

# ╔═╡ Cell order:
# ╠═49bb3630-4e0d-11ec-39f2-c3d1023de02a
# ╠═f1c7f04c-65e5-4cae-b3fe-42e3fefd058d
# ╠═d1651f70-47c1-4ce7-9f6c-bac9693d12aa
# ╠═e51ad0a1-60b6-4ac0-8a36-853cdf1b036a
# ╠═2a810db4-e771-432d-8880-45667807d507
# ╠═811f6a70-a9d6-41bb-bf39-d2a13c26b4cc
# ╠═5f8a9c94-b0b1-4bb3-ae81-a664863dbe53
# ╠═1e6a4ab2-837c-4aee-a432-0fec601b126f
# ╠═b1d4e2b4-d9ad-4ec7-91f7-39080b8bdd4b
# ╠═fe4b4a43-1c6f-4b6e-b21e-0938b41fa657
# ╠═4155d55a-ea9c-4343-ae56-69705b94984b
# ╠═1bc2e2c3-f7ef-44e8-832c-b239866f3c11
# ╠═e348d5ba-29a1-4f53-9608-d47d0afad7df
# ╠═ba00a0fa-a7de-493e-a86d-67d132281693
