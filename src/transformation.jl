"C-G coefficient c1 * c2 * <J, I, MJ, MI|J, I, F, MF>"
function *(b::Bra{T1,UncoupledHyperfineStructureState{L,S,J,I}}, k::Ket{T2,HyperfineStructureState{L,S,J,I}}) where {T1,T2,L,S,J,I}
    MJ, MI = b.s.MJ, b.s.MI
    F, MF = k.s.F, k.s.MF
    cg = clebschgordan(J, MJ, I, MI, F, MF)
    return b.c * k.c * cg
end

function *(b::Bra{T1,HyperfineStructureState{L,S,J,I}}, k::Ket{T2,UncoupledHyperfineStructureState{L,S,J,I}}) where {T1,T2,L,S,J,I}
    c = *(k', b')
    return c'
end

function *(bv::BraVec, kv::KetVec)
    out = 0
    for b in bv, k in kv
        out += b * k
    end
    return out
end

function basistransform(kv::KetVec{T,S1}, basis::AbstractVector{S2}) where {T,S1,S2}
    length(kv) == length(basis) || error("dimension not match")
    dims = length(basis)
    c = Vector{float(T)}(undef, dims)
    for i = 1:dims
        b = Bra(basis[i])
        cc = sum([b * k for k in kv])
        c[i] = abs(cc) > eps(float(T)) ? cc : zero(cc)
    end
    KetVec{float(T),S2}(c, basis)
end