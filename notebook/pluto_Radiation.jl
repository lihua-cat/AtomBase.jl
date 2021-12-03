### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ c12f4226-e920-48e8-9384-c5cec1413c2d
begin
	import Pkg
	Pkg.activate("..")
	Text(sprint(io->Pkg.status(io=io)))
end

# ╔═╡ 34451bc7-fed7-46a9-8609-c8e2c2a87a84
push!(LOAD_PATH, "C:\\Users\\lihua\\.julia\\dev")

# ╔═╡ ce85557a-79ae-4d28-8aaf-a3e1e57964c6
using Revise

# ╔═╡ 3950837f-4ab1-4ec8-beed-8092635517f1
using PlutoUI

# ╔═╡ d2bbd0b8-3862-403c-be88-48d2fc8cef08
using AtomBase

# ╔═╡ d64730eb-68c7-4f82-8a2a-d6d42ad4c419
using Unitful

# ╔═╡ b36dd160-d454-4e46-8894-c88e651c9e97
html"""<style>
main {
    max-width: 900px;
    align-self: flex-start;
    margin-left: 50px;
}
"""

# ╔═╡ c8de153d-4572-47a7-944b-f70673aab66c
md"## 1. The algebra of irreducible tensor operators"

# ╔═╡ 4204c151-1f8b-4cd2-8382-538835231655
@doc wigner_eckart

# ╔═╡ f0e7eb88-8b9a-4355-8c2c-1066143e26d4
@doc uncoup_T1

# ╔═╡ 65e742bd-b8a4-4c13-8a9c-38339c22d36f
@doc uncoup_T2

# ╔═╡ b7560a48-fc46-4038-ac74-7ed0c9682678
@doc reducedME

# ╔═╡ 310f23d5-7f36-4779-8898-71327fe8e2ba
md"## 2. Radiative transition matrix elements(E1 and M1)"

# ╔═╡ c3be6cf1-de4a-44b7-aa89-97c386557153
@doc relative_transitionME

# ╔═╡ 888d58a9-9fe0-4896-80fc-92fc25e8b3d9
@doc reducedME_E1

# ╔═╡ 58ae7db1-bd53-407d-9edf-53ab83e39de1
@doc reducedME_M1

# ╔═╡ fa73f020-9cfb-4d24-bd44-9fc863e6cb70
md"## 3. Einstein coefficients A and B"

# ╔═╡ 800a334e-6977-4439-bce7-e6c7088ffc6a
@doc aᵢⱼ

# ╔═╡ d6d5ba4e-9789-4c3e-8335-69cb5e3c901b
@doc einsteinA

# ╔═╡ f0a7642f-a8a5-44b3-aa88-76fb09737231
@doc σᵢⱼ

# ╔═╡ 1143243c-a2a1-4bb9-8ccf-395290547086
md"## Example: $^{127}\mathrm{I}$"

# ╔═╡ 4764c056-88d2-4223-aaf4-042fb8568c6f
begin
	L = 1
	S = 1/2
	J = (3/2, 1/2)
	I = 5/2
	F = (4, 3)
	MF = (3, 3)
end

# ╔═╡ 8154092d-239d-4cca-9f94-9229836205b2
rme = reducedME_M1(L, S, J[1], L, S, J[2])

# ╔═╡ cffd34a0-2505-4de1-988f-60bdabb19032
c = relative_transitionME(J[1], I, F[1], MF[1], J[2], I, F[2], MF[2], 1)

# ╔═╡ c2cec325-ca77-4d58-aecf-c1133455eb8d
s = (c * rme)^2

# ╔═╡ da54a4ff-a580-40d3-b1f6-98eefebc2add
k = 7603u"cm^-1"

# ╔═╡ 4e8a6397-a7d2-44a3-a58e-1a9b4836b754
a = aᵢⱼ(k, s)

# ╔═╡ 346a7803-d43b-4826-b49b-b53f3cfcdfa3
einsteinA(k, L, S, J[1], I, F[1], L, S, J[2], I, F[2], "M1")

# ╔═╡ Cell order:
# ╠═b36dd160-d454-4e46-8894-c88e651c9e97
# ╠═ce85557a-79ae-4d28-8aaf-a3e1e57964c6
# ╠═34451bc7-fed7-46a9-8609-c8e2c2a87a84
# ╠═c12f4226-e920-48e8-9384-c5cec1413c2d
# ╠═3950837f-4ab1-4ec8-beed-8092635517f1
# ╠═d2bbd0b8-3862-403c-be88-48d2fc8cef08
# ╠═d64730eb-68c7-4f82-8a2a-d6d42ad4c419
# ╟─c8de153d-4572-47a7-944b-f70673aab66c
# ╟─4204c151-1f8b-4cd2-8382-538835231655
# ╟─f0e7eb88-8b9a-4355-8c2c-1066143e26d4
# ╟─65e742bd-b8a4-4c13-8a9c-38339c22d36f
# ╟─b7560a48-fc46-4038-ac74-7ed0c9682678
# ╟─310f23d5-7f36-4779-8898-71327fe8e2ba
# ╟─c3be6cf1-de4a-44b7-aa89-97c386557153
# ╟─888d58a9-9fe0-4896-80fc-92fc25e8b3d9
# ╟─58ae7db1-bd53-407d-9edf-53ab83e39de1
# ╟─fa73f020-9cfb-4d24-bd44-9fc863e6cb70
# ╟─800a334e-6977-4439-bce7-e6c7088ffc6a
# ╟─d6d5ba4e-9789-4c3e-8335-69cb5e3c901b
# ╟─f0a7642f-a8a5-44b3-aa88-76fb09737231
# ╟─1143243c-a2a1-4bb9-8ccf-395290547086
# ╠═4764c056-88d2-4223-aaf4-042fb8568c6f
# ╠═8154092d-239d-4cca-9f94-9229836205b2
# ╠═cffd34a0-2505-4de1-988f-60bdabb19032
# ╠═c2cec325-ca77-4d58-aecf-c1133455eb8d
# ╠═da54a4ff-a580-40d3-b1f6-98eefebc2add
# ╠═4e8a6397-a7d2-44a3-a58e-1a9b4836b754
# ╠═346a7803-d43b-4826-b49b-b53f3cfcdfa3
