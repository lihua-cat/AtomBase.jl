function diagonal(op::Operator{T, S}) where {T, S}
    # ishermitian(op) || error("Operator maxtrix is not hermitian")
    (op.ks == op.bs && op.c ≈ op.c') || error("Operator maxtrix is not hermitian")
    m = op.c
    h = Hermitian(m)
    J = eltype(op.ks).parameters[3]
    I = eltype(op.ks).parameters[4]
    sortby = J >= I ? x -> x : x -> -x
    values, vectors = eigen(h, sortby = sortby)
    ket_vectors = Vector{KetVec{float(T), S}}(undef, size(op, 2))
    for i in 1:size(op, 1)
        vec = vectors[1, i] >= 0 ? vectors[:, i] : -vectors[:, i]
        ket_vectors[i] = KetVec(vec, op.ks)
    end
    return (;splits = values, states = ket_vectors)
end