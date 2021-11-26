## Angular Momentum Operators

# Jz
op_z(k::Ket, M::Symbol) = Ket(getfield(k.s, M) * k.c, k.s)

# J+, J-
function op_ladder(k::Ket{T, UncoupledHyperfineStructureState{L,S,J,I}}, M::Symbol, s::Int) where {T,L,S,J,I}
    m = getfield(k.s, M)
    if M == :MJ
        j = J
    elseif M == :MI
        j = I
    end
    abs(m) > j && error("|m| > j occurs")
    s in (+1, -1) || error("(s: $s) is not +1 or -1")
    j == s * m && return RationalRoot(0)
    c = k.c * RationalRoot(√((j - s * m) * (j + s * m + 1)))
    state = reinstantiate(k.s, M => m + s)
    return Ket(c, state)
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
    end
end

# composite operator
function ops(basis::AbstractVector{<:AtomState}, ops::Tuple)
    d = length(basis)
    ks = [Ket(s) for s in basis]
    c = ks' .* ∘(ops...).(ks)
    c = eltype(c) == Real ? float.(c) : c
    return Operator(Matrix(c'), basis, basis)
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