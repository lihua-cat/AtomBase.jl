function diagnoal(op::Operator{T, S}) where {T, S}
    ishermitian(op) || error("Operator maxtrix is not hermitian")
    m = op.c
    h = Hermitian(m)
    J = eltype(op.ks).parameters[3]
    I = eltype(op.ks).parameters[4]
    sortby = J >= I ? x -> x : x -> -x
    values, vectors = eigen(h, sortby = sortby)
    ket_vectors = Vector{KetVec{float(T), S}}(undef, size(op, 2))
    for i in 1:size(op, 1)
        ket_vectors[i] = KetVec(vectors[:, i], op.ks)
    end
    return (;values, vectors = ket_vectors)
end