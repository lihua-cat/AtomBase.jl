### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 92ecd7f0-46ab-11ec-3748-9da2e4977a67
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add(url = "https://github.com/lihua-cat/UsefulFunctions.jl")
	Pkg.add(url = "https://github.com/lihua-cat/AtomBase.jl")
	Pkg.add("DataFrames")
end

# ╔═╡ c37f5f93-3707-483b-b580-e5b5a189669f
using AtomBase

# ╔═╡ dd2ebed3-21b3-414e-b320-4696dfa75657
using UsefulFunctions 

# ╔═╡ becf314a-93a3-4e2e-a2c5-46755364130e
using DataFrames

# ╔═╡ d4872ed6-b4e5-4b20-a285-dd333aa3c0cc
md"# AtomBase"

# ╔═╡ 1fa72bdc-7bf0-4c93-a3cf-83df40ebb544
md"## Atom State"

# ╔═╡ a2f39c76-dbea-4662-a963-8fb844ce6c03
let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	F = 3
	MF = 2
	HyperfineStructureState{L,S,J,I}(F, MF)
end

# ╔═╡ 3fda0649-b1a2-46bd-8baf-bc9391d8bc9a
let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	MJ = -3//2
	MI = 3//2
	UncoupledHyperfineStructureState{L,S,J,I}(MJ, MI)
end

# ╔═╡ 7d2159ec-f201-4a0d-a8c8-8abe34bdff7e
md"## Dirac Notation"

# ╔═╡ fab4eacb-497c-4ab6-a2ee-baa77d218fda
let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	F = 3
	MF = 2
	s = HyperfineStructureState{L,S,J,I}(F, MF)
	k = Ket(s)
end

# ╔═╡ e78ea820-98c7-4c83-9e6e-62355ab582e8
let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	F = 3
	MF = 2
	s = HyperfineStructureState{L,S,J,I}(F, MF)
	Ket(s)'
end

# ╔═╡ 4dc97325-2e7b-4d35-b4ff-fcd76ec3e105
let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	MJ = -3//2
	MI = 3//2
	s = UncoupledHyperfineStructureState{L,S,J,I}(MJ, MI)
	Ket(s)
end

# ╔═╡ fc030395-109a-46ce-8e6f-9cd785bc6881
let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	MJ = -3//2
	MI = 3//2
	s = UncoupledHyperfineStructureState{L,S,J,I}(MJ, MI)
	Ket(s)'
end

# ╔═╡ 70cec9c6-c891-4ada-b354-a6ca5d6f7a73
md"**normalized**"

# ╔═╡ b29c0833-d5e7-402a-a086-2c5d5d1f850e
let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	F = 3
	MF = 2
	s = HyperfineStructureState{L,S,J,I}(F, MF)
	k = Ket(s)
	b = k'
	Text("$b * $k = $(b * k)")
end

# ╔═╡ a5228079-31bc-4307-a6ef-a95dea3384a8
md"**orthogonal**"

# ╔═╡ 5af92c73-cc1b-49b1-b782-b22e9e38e443
let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	F = 3
	MF = 2
	s = HyperfineStructureState{L,S,J,I}(F, MF)
	k = Ket(s)
	b = Ket(HyperfineStructureState{L,S,J,I}(F, MF-1))'
	Text("$b * $k = $(b * k)")
end

# ╔═╡ 8adc89f8-0b96-4cbe-932f-df1d43b115e6
md"**CG coefficient**"

# ╔═╡ 91c8844c-3265-4106-acb6-7eb19d206c01
let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	F = 3
	MF = 2
	s1 = HyperfineStructureState{L,S,J,I}(F, MF)
	k1 = Ket(s1)
	MJ = 1//2
	MI = 3//2
	s2 = UncoupledHyperfineStructureState{L,S,J,I}(MJ, MI)
	b2 = Ket(s2)'
	Text("$b2 * $k1 = $(b2 * k1)")
end

# ╔═╡ 4a9efa12-c06a-40f9-850d-f32d81c42a90
md"## Basis"

# ╔═╡ b688d23e-245d-46f5-b6b3-dfda4ea093e4
hfs_c = let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	basis_hfs(L, S, J, I, couple = true)
end

# ╔═╡ 1f99fbbe-22dd-4bd3-9e67-4697c9203233
hfs_uc = let
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
	basis_hfs(L, S, J, I, couple = false)
end

# ╔═╡ a1acd1ca-9918-444f-8343-4f7836280533
MF = 1

# ╔═╡ 1334a6d2-1e64-4ec2-a001-77f5495cd65e
basis_uc = hfs_uc[hfs_uc.MF .== MF, :basis1][]

# ╔═╡ 6e3ac2f9-3d1a-48c2-a536-7f0fed684e4a
basis_c = hfs_uc[hfs_uc.MF .== MF, :basis2][]

# ╔═╡ 887ae902-bd49-4270-a667-10a4409c28c0
md"## Operator Matrix Representation"

# ╔═╡ becc8535-a58d-452d-81c0-f5a17893d5f6
Base.float(op::Operator) = Operator(float.(op.c), op.ks, op.bs)

# ╔═╡ 59ad7b23-a03f-428d-af1a-1622ad44a0a5
Jz = AtomBase.𝐉𝐳(basis_uc)

# ╔═╡ e4bdef49-df82-4a08-b46c-95e4c8a4784f
J₊I₋ = AtomBase.𝐉₊𝐈₋(basis_uc)

# ╔═╡ 2a2146de-d99c-47ef-bfd0-adb4eb9977fe
J₋I₊ = AtomBase.𝐉₋𝐈₊(basis_uc)

# ╔═╡ 16444fe1-68a9-4b3f-b6e7-ef4a0eec514e
J₊²I₋² = AtomBase.𝐉₊²𝐈₋²(basis_uc)

# ╔═╡ 211d76e5-5e44-4cff-8191-27259bc7bfa8
md"## Diagnoalization"

# ╔═╡ 459f3741-3560-44e3-9c75-69f83a6293d7
let
	vals, vecs = diagonal(Jz)
end

# ╔═╡ fddba77b-cde9-4307-8167-fcd8833f7a5b
h = hamiltonian_hfs(basis_uc, 0.02759, 0.03812)

# ╔═╡ 919d240d-8e96-4235-a00e-c0d95d6a0d8c
vals, vecs = diagonal(h)

# ╔═╡ 55189e1d-0431-4317-ad4f-4ed71f679f4a
vecs2 = [basistransform(v, basis_c) for v in vecs]

# ╔═╡ ff6801e0-f10f-4dc0-8665-8cbef6cd28b5
kv = vecs[1]

# ╔═╡ 9f07f248-6fc1-433d-9358-30343feb4e54
md"## Basis Transformation"

# ╔═╡ 20f1d19d-7775-40b1-8952-20ca1d2137b9
kvt = basistransform(kv, basis_c)

# ╔═╡ 7b93d187-8ffa-41e2-9a48-506c997f3601
kvtt = basistransform(kvt, basis_uc)

# ╔═╡ 34e24f64-8bed-4268-ba0d-1c1c9fc5b205
kv.c ≈ kvtt.c

# ╔═╡ 8aaf45b8-a941-42c9-a8eb-f50d0e471496
md"## Radiation between HFS"

# ╔═╡ fe8a2c7d-86ed-4631-bf0b-8d8b0eedef8d
md"alkali D1 D2"

# ╔═╡ 432efe94-02de-4cfa-82a9-b2526aedec95
let
	L = (1, 0)
	S = (1//2, 1//2)
	J = (1//2, 1//2)
	I_list = 3//2:7//2
	transition = collect('a':'d')
	df1 = DataFrame(Transition = transition)
	for I in I_list
		p = []
		for Fl in I-J[2]:I+J[2]
			for Fu in I-J[1]:I+J[1]
				c = uncoup_T1(J[1], I, Fu, J[2], I, Fl, 1)^2
				push!(p, c)
			end
		end
		df1[!, Symbol("I=$I")] = p / sum(p)
	end
	df1
end

# ╔═╡ 8976e17c-f4f0-4e63-8e4f-8dc3820d812a
let
	L = (1, 0)
	S = (1//2, 1//2)
	J = (3//2, 1//2)
	I_list = 3//2:7//2
	transition = collect('e':'j')
	df2 = DataFrame(Transition = transition)
	for I in I_list
		p = []
		for Fl in I-J[2]:I+J[2]
			for Fu in Fl-1:Fl+1
				push!(p, uncoup_T1(J[1], I, Fu, J[2], I, Fl, 1)^2)
			end
		end
		df2[!, Symbol("I=$I")] = p / sum(p)
	end
	df2
end

# ╔═╡ 5def5d5c-c3af-4570-9b33-418149a72159
md"I127"

# ╔═╡ 0b5e2e9c-9058-447a-a983-e22c77ab71fc
begin
	L = (1, 1)
	S = (1//2, 1//2)
	J = (3//2, 1//2)
	I = 5//2
	Fl = I-J[1]:I+J[1]
	Fu = I-J[2]:I+J[2]
	df = DataFrame()
	F2v = []
	F1v = []
	p = []
	for F1 in Fl, F2 in Fu
		push!(F1v, Int(F1))
		push!(F2v, Int(F2))
		push!(p, uncoup_T1(J[1], I, F1, J[2], I, F2, 1)^2)
	end
	df.Fu = F2v
	df.Fl = F1v
	df.Relative = p / sum(p)
	df
end

# ╔═╡ 7d91c377-a5a3-43e6-babc-7fa7e871a11b
df[df.Fu .== 3 .&& df.Fl .== 4, :].Relative[1]

# ╔═╡ e5e05c16-c38a-4357-8376-bf9c1306d7b6
begin
	j = 3/2
	i = 1
	f = 5/2
	f*(f+1) - j*(j+1) - i*(i+1)
end

# ╔═╡ Cell order:
# ╟─d4872ed6-b4e5-4b20-a285-dd333aa3c0cc
# ╠═92ecd7f0-46ab-11ec-3748-9da2e4977a67
# ╠═c37f5f93-3707-483b-b580-e5b5a189669f
# ╠═dd2ebed3-21b3-414e-b320-4696dfa75657
# ╠═becf314a-93a3-4e2e-a2c5-46755364130e
# ╟─1fa72bdc-7bf0-4c93-a3cf-83df40ebb544
# ╟─a2f39c76-dbea-4662-a963-8fb844ce6c03
# ╟─3fda0649-b1a2-46bd-8baf-bc9391d8bc9a
# ╟─7d2159ec-f201-4a0d-a8c8-8abe34bdff7e
# ╟─fab4eacb-497c-4ab6-a2ee-baa77d218fda
# ╟─e78ea820-98c7-4c83-9e6e-62355ab582e8
# ╟─4dc97325-2e7b-4d35-b4ff-fcd76ec3e105
# ╟─fc030395-109a-46ce-8e6f-9cd785bc6881
# ╟─70cec9c6-c891-4ada-b354-a6ca5d6f7a73
# ╟─b29c0833-d5e7-402a-a086-2c5d5d1f850e
# ╟─a5228079-31bc-4307-a6ef-a95dea3384a8
# ╟─5af92c73-cc1b-49b1-b782-b22e9e38e443
# ╟─8adc89f8-0b96-4cbe-932f-df1d43b115e6
# ╟─91c8844c-3265-4106-acb6-7eb19d206c01
# ╟─4a9efa12-c06a-40f9-850d-f32d81c42a90
# ╟─b688d23e-245d-46f5-b6b3-dfda4ea093e4
# ╟─1f99fbbe-22dd-4bd3-9e67-4697c9203233
# ╠═a1acd1ca-9918-444f-8343-4f7836280533
# ╠═1334a6d2-1e64-4ec2-a001-77f5495cd65e
# ╠═6e3ac2f9-3d1a-48c2-a536-7f0fed684e4a
# ╟─887ae902-bd49-4270-a667-10a4409c28c0
# ╠═becc8535-a58d-452d-81c0-f5a17893d5f6
# ╠═59ad7b23-a03f-428d-af1a-1622ad44a0a5
# ╠═e4bdef49-df82-4a08-b46c-95e4c8a4784f
# ╠═2a2146de-d99c-47ef-bfd0-adb4eb9977fe
# ╠═16444fe1-68a9-4b3f-b6e7-ef4a0eec514e
# ╟─211d76e5-5e44-4cff-8191-27259bc7bfa8
# ╠═459f3741-3560-44e3-9c75-69f83a6293d7
# ╠═fddba77b-cde9-4307-8167-fcd8833f7a5b
# ╠═919d240d-8e96-4235-a00e-c0d95d6a0d8c
# ╠═55189e1d-0431-4317-ad4f-4ed71f679f4a
# ╠═ff6801e0-f10f-4dc0-8665-8cbef6cd28b5
# ╟─9f07f248-6fc1-433d-9358-30343feb4e54
# ╠═20f1d19d-7775-40b1-8952-20ca1d2137b9
# ╠═7b93d187-8ffa-41e2-9a48-506c997f3601
# ╠═34e24f64-8bed-4268-ba0d-1c1c9fc5b205
# ╟─8aaf45b8-a941-42c9-a8eb-f50d0e471496
# ╟─fe8a2c7d-86ed-4631-bf0b-8d8b0eedef8d
# ╟─432efe94-02de-4cfa-82a9-b2526aedec95
# ╟─8976e17c-f4f0-4e63-8e4f-8dc3820d812a
# ╟─5def5d5c-c3af-4570-9b33-418149a72159
# ╠═0b5e2e9c-9058-447a-a983-e22c77ab71fc
# ╠═7d91c377-a5a3-43e6-babc-7fa7e871a11b
# ╠═e5e05c16-c38a-4357-8376-bf9c1306d7b6
