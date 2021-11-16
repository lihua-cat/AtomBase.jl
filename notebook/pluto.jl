### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 92ecd7f0-46ab-11ec-3748-9da2e4977a67
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add(url = "https://github.com/lihua-cat/AtomBase.jl")
end

# ╔═╡ c37f5f93-3707-483b-b580-e5b5a189669f
using AtomBase

# ╔═╡ b688d23e-245d-46f5-b6b3-dfda4ea093e4
begin
	L = 1
	S = 1//2
	J = 3//2
	I = 5//2
end

# ╔═╡ 4c9f8926-f023-4d7d-ad8c-689143d57b8e
hfs_c = basis_hfs(L, S, J, I, couple = true)

# ╔═╡ 1f99fbbe-22dd-4bd3-9e67-4697c9203233
hfs_uc = basis_hfs(L, S, J, I, couple = false)

# ╔═╡ 1334a6d2-1e64-4ec2-a001-77f5495cd65e
basis_uc = basis_get(hfs_uc, :B1, :MF, 1)

# ╔═╡ 6e3ac2f9-3d1a-48c2-a536-7f0fed684e4a
basis_c = basis_get(hfs_uc, :B2, :MF, 1)

# ╔═╡ 59ad7b23-a03f-428d-af1a-1622ad44a0a5
Jz = 𝒥𝓏(basis_uc)

# ╔═╡ d8a5193e-b186-436d-a217-209b607a147d
Iz = ℐ𝓏(basis_uc)

# ╔═╡ e4bdef49-df82-4a08-b46c-95e4c8a4784f
J₊I₋ = 𝒥₊ℐ₋(basis_uc)

# ╔═╡ 312840f2-96f7-41d0-9ab8-edf2a8f8dae4
J₋I₊ = 𝒥₋ℐ₊(basis_uc)

# ╔═╡ 16444fe1-68a9-4b3f-b6e7-ef4a0eec514e
J₊²I₋² = 𝒥₊²ℐ₋²(basis_uc)

# ╔═╡ 459f3741-3560-44e3-9c75-69f83a6293d7
vals, vecs = diagonal(Jz)

# ╔═╡ f8ab3973-d277-47f4-90f3-0bbe34bb2ae9
kv = vecs[1]

# ╔═╡ 20f1d19d-7775-40b1-8952-20ca1d2137b9
kvt = basistransform(kv, basis_c)

# ╔═╡ 7b93d187-8ffa-41e2-9a48-506c997f3601
kvtt = basistransform(kvt, basis_uc)

# ╔═╡ Cell order:
# ╠═92ecd7f0-46ab-11ec-3748-9da2e4977a67
# ╠═c37f5f93-3707-483b-b580-e5b5a189669f
# ╠═b688d23e-245d-46f5-b6b3-dfda4ea093e4
# ╠═4c9f8926-f023-4d7d-ad8c-689143d57b8e
# ╠═1f99fbbe-22dd-4bd3-9e67-4697c9203233
# ╠═1334a6d2-1e64-4ec2-a001-77f5495cd65e
# ╠═6e3ac2f9-3d1a-48c2-a536-7f0fed684e4a
# ╠═59ad7b23-a03f-428d-af1a-1622ad44a0a5
# ╠═d8a5193e-b186-436d-a217-209b607a147d
# ╠═e4bdef49-df82-4a08-b46c-95e4c8a4784f
# ╠═312840f2-96f7-41d0-9ab8-edf2a8f8dae4
# ╠═16444fe1-68a9-4b3f-b6e7-ef4a0eec514e
# ╠═459f3741-3560-44e3-9c75-69f83a6293d7
# ╠═f8ab3973-d277-47f4-90f3-0bbe34bb2ae9
# ╠═20f1d19d-7775-40b1-8952-20ca1d2137b9
# ╠═7b93d187-8ffa-41e2-9a48-506c997f3601
