## Angular Momentum Operators

# Jz
op_z(k::Ket, M::Symbol) = Ket(getfield(k.s, M) * k.c, k.s)

# J+, J-
function op_ladder(k::Ket, M::Symbol, s::Int)
    m = getfield(k.s, M)
    if M == :MJ
        j = typeof(k.s).parameters[3]
    elseif M == :MI
        j = typeof(k.s).parameters[4]
    end
    abs(m) > j && return error("|m| > j occurs")
    if s == 1
        j == m && return zero(RationalRoot{Int})
        c = k.c * RationalRoot(√((j - m) * (j + m + 1)))
    elseif s == -1
        j == -m && return zero(RationalRoot{Int})
        c = k.c * RationalRoot(√((j + m) * (j - m + 1)))
    end
    s = reinstantiate(k.s, M => m + s)
    return Ket(c, s)
end

𝐉𝐳(k::Ket) = op_z(k, :MJ)
𝐉₊(k::Ket) = op_ladder(k, :MJ, +1)
𝐉₋(k::Ket) = op_ladder(k, :MJ, -1)
𝐈𝐳(k::Ket) = op_z(k, :MI)
𝐈₊(k::Ket) = op_ladder(k, :MI, +1)
𝐈₋(k::Ket) = op_ladder(k, :MI, -1)

# special case A⋅0 = 0 
for op in (:(𝐉𝐳), :(𝐉₊), :(𝐉₋), :(𝐈𝐳), :(𝐈₊), :(𝐈₋))
    @eval $(op)(x::Number) = x == 0 ? x : error("unexpected x: $x")
end

# single operator
for op in (:(𝐉𝐳), :(𝐉₊), :(𝐉₋), :(𝐈𝐳), :(𝐈₊), :(𝐈₋))
    @eval function $(op)(basis::AbstractVector{<:AtomState})
        d = length(basis)
        ks = [Ket(s) for s in basis]
        c = ks' .* $(op).(ks)
        return Operator(Matrix(c'), basis, basis)
        # return Operator(float.(c'), basis, basis)
    end
end

# composite operator
function ops(basis::AbstractVector{<:AtomState}, ops::Tuple)
    d = length(basis)
    ks = [Ket(s) for s in basis]
    c = ks' .* ∘(ops...).(ks)
    c = eltype(c) == Real ? float.(c) : c
    return Operator(Matrix(c'), basis, basis)
    # return Operator(float.(c'), basis, basis)
end

𝐉₊𝐈₋(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉₊, 𝐈₋))
𝐉₋𝐈₊(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉₋, 𝐈₊))

𝐉𝐳𝐉₊𝐈𝐳𝐈₋(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉𝐳, 𝐉₊, 𝐈𝐳, 𝐈₋))
𝐉𝐳𝐉₊𝐈₋𝐈𝐳(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉𝐳, 𝐉₊, 𝐈₋, 𝐈𝐳))
𝐉₊𝐉𝐳𝐈𝐳𝐈₋(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉₊, 𝐉𝐳, 𝐈𝐳, 𝐈₋))
𝐉₊𝐉𝐳𝐈₋𝐈𝐳(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉₊, 𝐉𝐳, 𝐈₋, 𝐈𝐳))

𝐉𝐳𝐉₋𝐈𝐳𝐈₊(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉𝐳, 𝐉₋, 𝐈𝐳, 𝐈₊))
𝐉𝐳𝐉₋𝐈₊𝐈𝐳(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉𝐳, 𝐉₋, 𝐈₊, 𝐈𝐳))
𝐉₋𝐉𝐳𝐈𝐳𝐈₊(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉₋, 𝐉𝐳, 𝐈𝐳, 𝐈₊))
𝐉₋𝐉𝐳𝐈₊𝐈𝐳(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉₋, 𝐉𝐳, 𝐈₊, 𝐈𝐳))

𝐉₊²𝐈₋²(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉₊, 𝐉₊, 𝐈₋, 𝐈₋))
𝐉₋²𝐈₊²(basis::AbstractVector{<:AtomState}) = ops(basis, (𝐉₋, 𝐉₋, 𝐈₊, 𝐈₊))