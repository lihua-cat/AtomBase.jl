using AtomBase
using HalfIntegers

L = 1
S = 1//2
J = 3//2
I = 5//2
hfs_c = basis_hfs(L, S, J, I, couple = true)
hfs_uc = basis_hfs(L, S, J, I, couple = false)

s1 = HyperfineStructureState{L,S,J,I}(2, 1)
s2 = HyperfineStructureState{L,S,J,I}(2, 2)
s1 == s2
Ket(s1) + Ket(s2)
sv = [HyperfineStructureState{L,S,J,I}(2, i) for i = 2:-1:-2]
ke = Ket(1//2, s1)

MF = 1
basis = hfs_uc[hfs_uc.MF .== 1, :basis1][]
basis2 = hfs_uc[hfs_uc.MF .== 1, :basis2][]

Jz = AtomBase.𝐉𝐳(basis)
Iz = AtomBase.𝐈𝐳(basis)
J₊I₋ = AtomBase.𝐉₊𝐈₋(basis)
J₋I₊ = AtomBase.𝐉₋𝐈₊(basis)
J₊²I₋² = AtomBase.𝐉₊²𝐈₋²(basis)

vals, vecs = diagonal(Jz)
vals, vecs = diagonal(J₊I₋ + J₋I₊)

kv = vecs[1]
kvt = basistransform(kv, basis2)
kvtt = basistransform(kvt, basis)
kvtt.c ≈ kv.c

J = (3/2, 1/2)
I = (5/2, 5/2)
for F1 in J[1]+I[1]:-1:abs(J[1]-I[1]), F2 in J[2]+I[2]:-1:abs(J[2]-I[2])
    c = uncoup_T1(J[1], I[1], F1, J[2], I[2], F2, 1)^2
    println("$F2 -> $F1: $c")
end

h = hamiltonian_hfs(basis, 0.02759, 0.03812)
h.c
vals, vecs = diagonal(h)
vals
vecs2 = [basistransform(v, basis2) for v in vecs]

k = Ket(HyperfineStructureState{1,1/2,1/2,5/2}(3, 3))
b = Bra(HyperfineStructureState{1,1/2,3/2,5/2}(4, 4))

transitionME(3/2, 5/2, 4, 4, 1/2, 5/2, 3, 3, 1)

relative_transition_intensity(b, k)