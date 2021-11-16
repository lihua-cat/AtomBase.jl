using AtomBase

L = 1
S = 1//2
J = 3//2
I = 5//2
hfs_c = basis_hfs(L, S, J, I, couple = true)
hfs_uc = basis_hfs(L, S, J, I, couple = false)

s1 = HyperfineStructureState{L,S,J,I}(2, 1)
s2 = HyperfineStructureState{L,S,J,I}(2, 2)
s1 == s2
sv = [HyperfineStructureState{L,S,J,I}(2, i) for i = 2:-1:-2]
ke = Ket(1//2, s1)

basis = basis_get(hfs_uc, :B1, :MF, 1)
basis2 = basis_get(hfs_uc, :B2, :MF, 1)

Jz = 𝒥𝓏(basis)
Iz = ℐ𝓏(basis)
J₊I₋ = 𝒥₊ℐ₋(basis)
J₋I₊ = 𝒥₋ℐ₊(basis)
J₊²I₋² = 𝒥₊²ℐ₋²(basis)

vals, vecs = diagonal(Jz)
vals, vecs = diagonal(J₊I₋ + J₋I₊)

kv = vecs[1]
kvt = basistransform(kv, basis2)
kvtt = basistransform(kvt, basis)
kvtt.c ≈ kv.c

J = (1/2, 1/2)
I = (7/2, 7/2)
for F1 in J[1]+I[1]:-1:abs(J[1]-I[1]), F2 in J[2]+I[2]:-1:abs(J[2]-I[2])
    c = uncoup_T1(J[1], I[1], F1, J[2], I[2], F2, 1)^2
    println("$F1 -> $F2: $c")
end