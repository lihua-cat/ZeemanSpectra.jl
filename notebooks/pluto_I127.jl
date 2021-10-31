### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 975815a4-7825-4cb0-8628-c8c4378221e8
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add("DataFrames")
	Pkg.add("PhysicalConstants")
	Pkg.add("Unitful")
	Pkg.add("HalfIntegers")
	Pkg.add("UnitfulRecipes")
	Pkg.add("PlutoUI")
	Pkg.add("Plots")
	Pkg.add("LaTeXStrings")
	Pkg.add(url = "https://github.com/lihua-cat/LineProfile.jl")
	Pkg.add(url = "https://github.com/lihua-cat/ZeemanSpectra.jl")
end

# ╔═╡ b1bb7864-18c3-4074-bf78-90228b11721d
# import ZeemanSpectra: atom_df, zeeman_struc, zeeman_spec, profile_voigt, fwhm_doppler, transition_matrix_element_m1_total, k_I127, A_I127, σ0_I127
using ZeemanSpectra

# ╔═╡ 0ae40869-59a2-4e54-b205-61f52c6480d8
using HalfIntegers

# ╔═╡ f9b30b0c-f219-40f4-ac2c-6eefcb5f12d6
using DataFrames

# ╔═╡ be11d03a-0646-43d1-a1f6-186e6a5bdb87
using Unitful, UnitfulRecipes

# ╔═╡ b7a234a2-b520-4417-8814-5f5d51723521
using LaTeXStrings

# ╔═╡ a2d38136-b186-4984-a960-ced4effad6cf
using Plots;gr()

# ╔═╡ bc55edf6-87c7-4408-8808-a42ac364f916
using PlutoUI

# ╔═╡ 350d9f41-85f8-457c-8b5a-a6a86ce95d62
md"# Zeeman Spectra"

# ╔═╡ c3023a88-bbd7-4ae5-b3bf-1c799bd5cec0
import PhysicalConstants.CODATA2018: c_0, k_B, h, μ_0, μ_B

# ╔═╡ f7bc1fc7-9fbc-4189-a304-6e3e78b6af67
default(;framestyle=:box, widen=true, fg_legend = :gray80, dpi=300, gridalpha=0.25, fontfamily="Source Han Sans")

# ╔═╡ cf7b4b7b-e889-483c-98f0-3b30ce75b595
md"# 1. atom data"

# ╔═╡ bc4c88ac-d80f-4bcd-bc4d-341bbfb07c7c
ATOM_DATA

# ╔═╡ 1c015f1b-f0b3-4f76-836b-2d6f40f0792b
md"# 2. Zeeman splits"

# ╔═╡ e6a16e4d-91b8-483e-aa0f-77c2a4a5f6d9
atom_name = "I127"

# ╔═╡ dd73de36-5be8-4f8c-9f40-c81ca1c3bbf2
atom_state = "g"

# ╔═╡ b3d6ba46-2f02-44f2-a34b-a8038b0f882b
df0 = zeeman_struc(atom_name, atom_state)

# ╔═╡ 9ec96c32-efa0-458b-a405-5f7a25334a45
B_range = (0:1:3000)u"Gauss"

# ╔═╡ a42df08c-5b2e-4e28-a842-9860c36b2591
begin
	F0 = [4, 3, 2, 1]
	df00 = filter(row -> row.F in F0, df0)
	Vec_F = df00[!, :F]
	Vec_MF = df00[!, :MF]
	Mat_E = []
	for B in B_range
    	df01 = zeeman_struc(atom_name, atom_state, B)
		filter!(row -> row.F in F0, df01)
    	Mat_E = isempty(Mat_E) ? df01[!, :E] : hcat(Mat_E, df01[!, :E]) 
	end
end

# ╔═╡ 27387dc4-76a3-4cc5-8c90-e2c5d57c4119
begin
	plot()
	for nrow in 1 : size(Mat_E)[1]
		plot!(B_range, 
			Mat_E[nrow, :], 
			label = "F=$(Vec_F[nrow]) MF =$(Vec_MF[nrow])", lw = 1, 
			legend = :outerright, legendfontsize = 6)
	end
	title!("Zeeman split")
	xlabel!("Ambient magnetic field(Gauss)")
	ylabel!("Energy level(cm^-1)")
end

# ╔═╡ 112b2b7b-549f-424c-8bdd-e611306e2480
md"# 3. Transition cross section"

# ╔═╡ c26b2e20-f568-4e35-823f-4871dc296e46
atom_state1 = "e1"

# ╔═╡ 2d77fc21-8096-44a9-b500-7eddf64c499b
atom_state2 = "g"

# ╔═╡ 3df50cb9-95a7-43dd-a19f-1c4222a8f557
@bind T Slider(0:400, default = 200, show_value=true)

# ╔═╡ 15376045-e511-4011-b304-ed84ee537b0e
@bind P Slider(0:20, default = 10, show_value=true)

# ╔═╡ 21ee843a-9f6a-4193-a761-0449a7534a14
νp = 2 * 5.0 * P

# ╔═╡ 7de22e25-ddd4-41f9-a469-bd3321485e8b
md"### 1. no magentic field"

# ╔═╡ 5457ab5f-faa5-46e5-8479-f4332b4c0f8d
begin 
	df1 = zeeman_struc(atom_name, atom_state1)
	df2 = zeeman_struc(atom_name, atom_state2)
	k34_0 = abs(df1[1, :E] - df2[1, :E])
end

# ╔═╡ ac3e2486-b453-475f-a763-f83d3c6dd13c
νd = fwhm_doppler(k34_0*c_0, ATOM_DATA[ATOM_DATA.Name .== atom_name, :M][], T*u"K") |> u"MHz"

# ╔═╡ 6866a683-b30a-457f-86cd-d53c5c65f8f6
kx = collect(7602.2:0.0001:7603.8)u"cm^-1"

# ╔═╡ bb71c05f-151d-4078-91c6-27791ab46d5f
#einstein_A_m1(λ, L, S, J, I, F) = uconvert(u"s^-1", 64π^4 / (3h * λ^3) * μ_0 / 4π * μ_B^2 * transition_matrix_element_m1_total(L, S, J, I, F)^2) / (2F[1] + 1)

# ╔═╡ 61fbcc32-0d38-47df-9ae5-5a28e21b9d25
line_profile_peak = uconvert(u"s", profile_voigt(0.0u"MHz";ν0 = 0.0u"MHz", νd = νd, νp = νp*u"MHz"))

# ╔═╡ 53509b73-044d-4412-afe8-2e595a69bbfa
F1_list = 3:-1:2

# ╔═╡ 91a2ddee-d130-45e9-b12a-6b1deebd8721
F2_list = 4:-1:1

# ╔═╡ c0425020-25cc-4ce9-9a7a-4766824831b4
begin
	L = (1, 1)
	S = (1/2, 1/2)
	J = (1/2, 3/2)
	I = (5/2, 5/2)
	df_A = DataFrame(F1 = HalfInt[], F2 = HalfInt[], k = Unitful.Wavenumber[], A = Unitful.Frequency[], gA = Unitful.Frequency[], σm = Unitful.Area[], tme = Rational[], σ = Vector{<:Unitful.Area}[])
	for F1 in F1_list
		for F2 in F2_list
			k = abs(df1[df1.F .== F1, :E][1] - df2[df2.F .== F2, :E][1])
			A = einstein_A_m1(k, L, S, J, I, (F1, F2))
			gA = (2F1 + 1) * A
			σm = A / k^2 / 8π * line_profile_peak
			σ = A / k^2 / 8π * uconvert.(u"s", profile_voigt.(u"MHz".(kx.*c_0);ν0 = u"MHz"(k*c_0), νd = νd, νp = νp*u"MHz"))
			tme = transition_matrix_element_m1_total(L, S, J, I, (F1, F2))^2
			push!(df_A, (F1, F2, k, A, gA, σm, tme, σ))
		end
	end
end

# ╔═╡ ae2dabb6-bc67-4caa-8adc-db01f9c317cc
begin
	df_A.σr = df_A.gA / sum(df_A[!, :gA])
	df_A.relative = df_A.tme / sum(df_A[!, :tme])
	select!(df_A, Not(:tme))
	df_A
end

# ╔═╡ a1b1392a-6b27-42c3-8c46-9f9639bcaaf6
σ0_ls = σ0_I127(3, 4, T*u"K", P*u"Torr", 5u"MHz/Torr")[2]

# ╔═╡ 63e89889-143d-443c-a79b-4c4d3c280d36
let
	plot(legend = false)
	x = ustrip.(u"cm^-1", kx)
	color = [:red, :blue, :yellow, :green, :pink, :purple, :cyan, :orange]
	σt = zeros(length(kx))
	i = 1
	for row in eachrow(df_A)
		σ = ustrip.(NoUnits, row.σ/σ0_ls)
		k = ustrip.(u"cm^-1", row.k)
		σm = ustrip(NoUnits, row.σm/σ0_ls)
		plot!(x, σ, lw = 1, color = color[i])
		annotate!([(k, σm, Plots.text("$(row.F1)-$(row.F2)", 8, color[i], :bottom))])
		i += 1
		σt .+= σ
	end
	# plot!(x, σt, lw = 1, color = :gray, ls = :dash)
	title!(L"\textrm{\sffamily Zeeman spectra}")
	ylabel!(L"\sigma/\sigma_0")
	xlabel!(L"k(\mathrm{cm^{-1}})")
end

# ╔═╡ be9d0b08-ee75-4f55-bea4-363f2cb31eaf
md"### 2. with magentic field"

# ╔═╡ da670bf1-82a8-4575-8a30-e0be282291c0
@bind BF Slider(0:1000, default = 200, show_value=true)

# ╔═╡ c409ed85-e506-4740-93d9-68f9ac4fb060
kc = uconvert(u"cm^-1", νp*u"MHz"/c_0)

# ╔═╡ 3bd2fbac-8895-417e-a923-b64aac906f78
upreferred(νp*u"MHz"/νd)

# ╔═╡ fab9a0c5-3255-4dc1-af66-e7f912746f93
begin
	iodine_spec(kx, T, νp, BF) = zeeman_spec(atom_name, atom_state1, atom_state2, kx, T*u"K", νp*u"MHz", BF*u"Gauss")
	df_E, df_CS = iodine_spec(kx, T, νp, BF)
	df_CS
end

# ╔═╡ e0745cb9-5a23-4528-a2d6-9e8b4073f1d9
df_F1_F2(F1, F2, pl, df = df_CS) = filter(row -> row.polarization in pl, groupby(groupby(df, :F1)[(F1=F1,)], :F2)[(F2=F2,)])	

# ╔═╡ 1ceb5117-5d8d-4915-b9ae-67e456114fa7
df34 = df_F1_F2(3, 4, ("S",))

# ╔═╡ d9d9dc2e-95ae-4c4d-8f52-f8b5f59e9aa0
begin
	df_pl_S = DataFrame(F1 = HalfInt[], F2 = HalfInt[], km= Unitful.Wavenumber[], σ = Vector{<:Unitful.Area}[], σm = Unitful.Area[])
	for F1 in F1_list
		for F2 in F2_list
			σ = sum(df_F1_F2(F1, F2, ("S",))[!, :σ])/2*3
			σm, km_i = findmax(σ)
			push!(df_pl_S, (F1, F2, kx[km_i], σ, σm))
		end
	end
	df_pl_S
end

# ╔═╡ d4356517-fda4-4169-9a87-928b7cf416b3
let
	plot(legend = false)
	x = ustrip.(u"cm^-1", kx)
	color = [:red, :blue, :yellow, :green, :pink, :purple, :cyan, :orange]
	i = 1
	ymax= 0
	for row in eachrow(df_pl_S)
		σ = ustrip.(NoUnits, row.σ/σ0_ls/7)
		k = ustrip.(u"cm^-1", row.km)
		σm = ustrip(NoUnits, row.σm/σ0_ls/7)
		plot!(x, σ, lw = 1, color = color[i])
		annotate!([(k, σm, Plots.text("$(row.F1)-$(row.F2)", 8, color[i], :bottom))])
		i += 1
		ymax = ymax > σm ? ymax : σm
	end
	title!(L"\textrm{\sffamily Zeeman spectra (S)}")
	annotate!([(7602.25, 0.9ymax, Plots.text("B = $BF", 16, :blue, :left))])
	ylabel!(L"\sigma/\sigma_0")
	xlabel!(L"k(\mathrm{cm^{-1}})")
end

# ╔═╡ 089b712c-d1aa-4413-8dca-05d07039abcc
let
	plot(legend = false)
	kx = collect(7603.07:0.0001:7603.2)u"cm^-1"
	df_E, df_CS = iodine_spec(kx, T, νp, BF)
	df34 = df_F1_F2(3, 4, ("S",), df_CS)
	x = ustrip.(u"cm^-1", kx)
	i = 1
	σt = zeros(length(kx))
	for row in eachrow(df34)
		σ = ustrip.(NoUnits, row.σ/σ0_ls/7)
		k = ustrip.(u"cm^-1", row.k0)
		σm = ustrip(NoUnits, row.σ0*line_profile_peak/σ0_ls/7)
		plot!(x, σ, lw = 1)
		annotate!([(k, σm, Plots.text("$(row.MF1)->$(row.MF2)", 6, :bottom))])
		σt .+= σ
		i += 1
	end
	plot!(x, σt, lw = 1, color = :gray, ls = :dash)
	ymax = maximum(σt)
	title!(L"\textrm{\sffamily Zeeman spectra (S)}")
	annotate!([(7602.25, 0.9ymax, Plots.text("B = $BF", 16, :blue, :left))])
	ylabel!(L"\sigma/\sigma_0")
	xlabel!(L"k(\mathrm{cm^{-1}})")
end

# ╔═╡ 6f336189-f859-4714-906f-48446a46e7b3
function σ_B(B_range, F1_list, F2_list)
	σ_ma = Matrix{Unitful.Area}(undef, length(B_range), length(F1_list)*length(F2_list))
	for i in 1:length(B_range)
		df_E, df_CS = iodine_spec(kx, T, νp, B_range[i])
		j = 1
		for F1 in F1_list
			for F2 in F2_list
				σ = sum(df_F1_F2(F1, F2, ("S",), df_CS)[!, :σ])/2*3
				σm = maximum(σ)
				σ_ma[i, j] = σm
				j += 1
			end
		end
	end
	return σ_ma
end

# ╔═╡ b6104955-ece9-4b55-9572-f3c300d9b78b
begin
	B_range2 = 0:5:400
	σ_ma = σ_B(B_range2, F1_list, F2_list)
end

# ╔═╡ 77e24d0b-83d9-4567-8cf5-d075b5bf7b14
let
	plot()
	x = ustrip.(u"cm^-1", kx)
	color = [:red, :blue, :yellow, :green, :pink, :purple, :cyan, :orange]
	i = 1
	for F1 in F1_list, F2 in F2_list
		σm = ustrip.(NoUnits, σ_ma[:, i]/σ0_ls)
		plot!(B_range2, σm, lw = 1, color = color[i], label = "F'=$F1->F=$F2")
		i += 1
	end
	title!(L"\textrm{\sffamily Zeeman spectra (S)}")
	ylabel!(L"\sigma(10^{-16}\mathrm{cm^2})")
	xlabel!(L"B(\mathrm{Gauss})")
end

# ╔═╡ f4a173ae-71dc-407f-8420-ae59c2f6a22f
function plot_spec(BF; kx = kx, atom_name = atom_name, atom_state1 = atom_state1, 
	atom_state2 = atom_state2, T = T, νp = νp)
	df_E, df_CS = zeeman_spec(atom_name, atom_state1, atom_state2, kx, T*u"K", νp*u"MHz", BF*u"Gauss")
	plot()
	yp = sum(df_CS[df_CS.polarization .== "P", :σ])*1e18
	ys = sum(df_CS[df_CS.polarization .== "S", :σ])/2*1e18
	ymax = ustrip(unit(yp[1]), maximum(vcat(yp, ys)))
	plot!(kx, yp, lw = 1, label = "P")
	plot!(kx, ys, lw = 1, label = "S")
	annotate!([(7602.25, 0.9ymax, Plots.text("B = $BF", 16, :blue, :left))])
	annotate!([(7602.25, 0.8ymax, Plots.text("νp = $(round(ustrip(u"cm^-1", kc), sigdigits = 2))"*" cm^-1", 8, :red, :left))])
	annotate!([(7602.25, 0.75ymax, Plots.text("νd = $(round(ustrip(u"cm^-1", νd/c_0), sigdigits = 2))"*" cm^-1", 8, :red, :left))])
	title!(L"\textrm{\sffamily Zeeman spectra}")
	ylabel!(L"\textrm{\sffamily transition cross section}(10^{-18}\mathrm{cm^2})")
	xlabel!(L"\textrm{\sffamily Wavenumber}(\mathrm{cm^{-1}})")
end

# ╔═╡ d84030c2-c707-4129-94f3-4709bc2d82e5
plot_spec(BF);plot!(;legend = :topright)

# ╔═╡ Cell order:
# ╟─350d9f41-85f8-457c-8b5a-a6a86ce95d62
# ╠═975815a4-7825-4cb0-8628-c8c4378221e8
# ╠═b1bb7864-18c3-4074-bf78-90228b11721d
# ╠═c3023a88-bbd7-4ae5-b3bf-1c799bd5cec0
# ╠═0ae40869-59a2-4e54-b205-61f52c6480d8
# ╠═f9b30b0c-f219-40f4-ac2c-6eefcb5f12d6
# ╠═be11d03a-0646-43d1-a1f6-186e6a5bdb87
# ╠═b7a234a2-b520-4417-8814-5f5d51723521
# ╠═a2d38136-b186-4984-a960-ced4effad6cf
# ╠═f7bc1fc7-9fbc-4189-a304-6e3e78b6af67
# ╠═bc55edf6-87c7-4408-8808-a42ac364f916
# ╟─cf7b4b7b-e889-483c-98f0-3b30ce75b595
# ╠═bc4c88ac-d80f-4bcd-bc4d-341bbfb07c7c
# ╟─1c015f1b-f0b3-4f76-836b-2d6f40f0792b
# ╟─e6a16e4d-91b8-483e-aa0f-77c2a4a5f6d9
# ╟─dd73de36-5be8-4f8c-9f40-c81ca1c3bbf2
# ╠═b3d6ba46-2f02-44f2-a34b-a8038b0f882b
# ╟─9ec96c32-efa0-458b-a405-5f7a25334a45
# ╠═a42df08c-5b2e-4e28-a842-9860c36b2591
# ╟─27387dc4-76a3-4cc5-8c90-e2c5d57c4119
# ╟─112b2b7b-549f-424c-8bdd-e611306e2480
# ╠═c26b2e20-f568-4e35-823f-4871dc296e46
# ╠═2d77fc21-8096-44a9-b500-7eddf64c499b
# ╠═3df50cb9-95a7-43dd-a19f-1c4222a8f557
# ╠═15376045-e511-4011-b304-ed84ee537b0e
# ╠═21ee843a-9f6a-4193-a761-0449a7534a14
# ╠═ac3e2486-b453-475f-a763-f83d3c6dd13c
# ╟─7de22e25-ddd4-41f9-a469-bd3321485e8b
# ╠═5457ab5f-faa5-46e5-8479-f4332b4c0f8d
# ╠═6866a683-b30a-457f-86cd-d53c5c65f8f6
# ╠═bb71c05f-151d-4078-91c6-27791ab46d5f
# ╠═61fbcc32-0d38-47df-9ae5-5a28e21b9d25
# ╠═53509b73-044d-4412-afe8-2e595a69bbfa
# ╠═91a2ddee-d130-45e9-b12a-6b1deebd8721
# ╠═c0425020-25cc-4ce9-9a7a-4766824831b4
# ╠═ae2dabb6-bc67-4caa-8adc-db01f9c317cc
# ╠═a1b1392a-6b27-42c3-8c46-9f9639bcaaf6
# ╟─63e89889-143d-443c-a79b-4c4d3c280d36
# ╟─be9d0b08-ee75-4f55-bea4-363f2cb31eaf
# ╠═da670bf1-82a8-4575-8a30-e0be282291c0
# ╠═c409ed85-e506-4740-93d9-68f9ac4fb060
# ╠═3bd2fbac-8895-417e-a923-b64aac906f78
# ╠═fab9a0c5-3255-4dc1-af66-e7f912746f93
# ╠═e0745cb9-5a23-4528-a2d6-9e8b4073f1d9
# ╠═1ceb5117-5d8d-4915-b9ae-67e456114fa7
# ╠═d9d9dc2e-95ae-4c4d-8f52-f8b5f59e9aa0
# ╠═d4356517-fda4-4169-9a87-928b7cf416b3
# ╠═089b712c-d1aa-4413-8dca-05d07039abcc
# ╟─6f336189-f859-4714-906f-48446a46e7b3
# ╠═b6104955-ece9-4b55-9572-f3c300d9b78b
# ╟─77e24d0b-83d9-4567-8cf5-d075b5bf7b14
# ╠═f4a173ae-71dc-407f-8420-ae59c2f6a22f
# ╟─d84030c2-c707-4129-94f3-4709bc2d82e5
