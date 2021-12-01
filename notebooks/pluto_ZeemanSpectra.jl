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

# ╔═╡ 55c722e9-6cf0-465a-bd63-aa15d417ae34
using Printf

# ╔═╡ 2a810db4-e771-432d-8880-45667807d507
import PhysicalConstants.CODATA2018: h, μ_B, μ_0

# ╔═╡ 8d1d5ad5-5f2b-4ffc-8cc5-959aa79ba342
md"## Atom Data"

# ╔═╡ 5f8a9c94-b0b1-4bb3-ae81-a664863dbe53
atom = I127

# ╔═╡ a2bafc87-52ef-4667-b907-9c9d920dd258
md"## Fine Structure"

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
md"## Hyperfine structure"

# ╔═╡ 73461ba8-9a53-4201-8022-57fa1d29a173
md"### Energy level splits"

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
md"### Transition lines"

# ╔═╡ 48777912-33ef-474d-bcdc-13caac6bb977
md"1.wavenumber, relative intensity and Einstein A coefficient"

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
		A = 64π^4 / 3h * k^3 * c / (2F2 + 1) |> u"s^-1"
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
	``^2P_{1/2} \rightarrow ^2P_{3/2}``: $(A)

	A in [3]: 7.87 s^-1
	
	k and c: see Table 4 in [1]
	"""
end

# ╔═╡ df5f0239-6f1a-4e6b-9236-9738dbc20fdd
uncoup_T1(3/2, 5/2, 4, 1/2, 5/2, 3, 1)^2 * 7.939805013011782 * 2 / 7

# ╔═╡ 15417b8a-fd8c-4174-9bea-7aadd422660f
uncoup_T1(3/2, 5/2, 3, 1/2, 5/2, 3, 1)^2 * 7.939805013011782 * 2 / 7

# ╔═╡ 70d1790d-7bb9-456d-82a3-cadae892803c
md"2.transition cross section"

# ╔═╡ 4155d55a-ea9c-4343-ae56-69705b94984b
df_spec = let
	BF = 0.0u"Gauss"
	kx = collect(7602.2:0.0001:7603.8)u"cm^-1"
	
	T = 200u"K"
	P = 0u"Torr"
	νp = 2 * 5u"MHz/Torr" * P
	
	df_spec = zeeman_spec(atom, 0, 1, kx, T, νp, "M1", BF)
	df = DataFrame(F1=[], F2=[], k=[], c=[], A=[], σ0=[], σ=[])

	for F1 in Set(df_spec.F1), F2 in Set(df_spec.F2)
		df34 = filter(row->row.F1 == F1 && row.F2 == F2, df_spec)
		c = sum(df34.c0)
		c < eps() && continue
		@assert all(df34.k0 .≈ df34.k0[1]) "error"
		A = sum(df34.a)/(2F2+1)
		σ0 = sum(df34.σ0)/(2F2+1)
		σ = sum(df34.σ)/(2F2+1)
		push!(df, (F1, F2, df34.k0[1], c, A, σ0, σ))
	end
	sort!(df, [:c, :A], rev=true)
	df
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
σ0_I127(4, 3, T = 200u"K", P = 0u"Torr") * 7/12

# ╔═╡ f33dfb79-0b28-4d14-95b5-4c49c075114c
md"""
## Reference
1. Engleman, R., Keller, R. A. & Palmer, B. A. Hyperfine structure and isotope shift of the 1.3-μm transition of ^129l. Appl Optics 19, 2767 (1980).
2. Kelly, M. A., McIver, J. K., Shea, R. F. & Hager, G. D. Frequency tuning of a CW atomic iodine laser via the Zeeman effect. Ieee J Quantum Elect 27, 263–273 (1991).
3. Lilenfeld, H. V., Richardson, R. J. & Hovis, F. E. The electron spin resonance spectrum of the 2 P 1/2 state of atomic iodine. J Chem Phys 74, 2129–2132 (1981).
4. Hager, G. D. et al. A simplified analytic model for gain saturation and power extraction in the flowing chemical oxygen-iodine laser. Ieee J Quantum Elect 32, 1525–1536 (1996).
  
"""

# ╔═╡ Cell order:
# ╠═49bb3630-4e0d-11ec-39f2-c3d1023de02a
# ╠═f1c7f04c-65e5-4cae-b3fe-42e3fefd058d
# ╠═d1651f70-47c1-4ce7-9f6c-bac9693d12aa
# ╠═e51ad0a1-60b6-4ac0-8a36-853cdf1b036a
# ╠═2a810db4-e771-432d-8880-45667807d507
# ╠═811f6a70-a9d6-41bb-bf39-d2a13c26b4cc
# ╠═55c722e9-6cf0-465a-bd63-aa15d417ae34
# ╟─8d1d5ad5-5f2b-4ffc-8cc5-959aa79ba342
# ╟─5f8a9c94-b0b1-4bb3-ae81-a664863dbe53
# ╟─a2bafc87-52ef-4667-b907-9c9d920dd258
# ╟─61054b00-6b53-4fab-abd4-f454f80f59e1
# ╟─773bab65-0dcd-4400-9b7f-8dec6effad8e
# ╟─e0b645bd-7909-422e-8da1-b6c340926a8f
# ╟─73461ba8-9a53-4201-8022-57fa1d29a173
# ╟─3cfa4727-37df-41c8-95fb-bcbd62fa047e
# ╟─537cf404-25f3-49b5-a9b6-09f68384e275
# ╟─abbd971e-f498-4f1f-bc88-f6064c3eee64
# ╟─430de0d6-442d-49c1-a2d7-6b22b8b99afc
# ╟─8ade433b-9039-4ef6-9a40-5f219e03b93a
# ╟─2cbf8cb2-95c8-4ed3-b82f-f59e40693ab8
# ╟─48777912-33ef-474d-bcdc-13caac6bb977
# ╟─e348d5ba-29a1-4f53-9608-d47d0afad7df
# ╟─5b049fac-21f3-4053-8d09-da7ae5a8fd22
# ╠═df5f0239-6f1a-4e6b-9236-9738dbc20fdd
# ╠═15417b8a-fd8c-4174-9bea-7aadd422660f
# ╟─70d1790d-7bb9-456d-82a3-cadae892803c
# ╟─4155d55a-ea9c-4343-ae56-69705b94984b
# ╟─cd826b82-1249-4e88-ac2b-7a25bfe24f61
# ╠═5112acd1-614a-48d3-b826-9a6a56b3ade9
# ╟─f33dfb79-0b28-4d14-95b5-4c49c075114c
