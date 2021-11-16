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

𝒥𝓏(k::Ket) = op_z(k, :MJ)
𝒥₊(k::Ket) = op_ladder(k, :MJ, +1)
𝒥₋(k::Ket) = op_ladder(k, :MJ, -1)
ℐ𝓏(k::Ket) = op_z(k, :MI)
ℐ₊(k::Ket) = op_ladder(k, :MI, +1)
ℐ₋(k::Ket) = op_ladder(k, :MI, -1)

# special case A⋅0 = 0 
for op in (:(𝒥𝓏), :(𝒥₊), :(𝒥₋), :(ℐ𝓏), :(ℐ₊), :(ℐ₋))
    @eval $(op)(x::Number) = x == 0 ? x : error("unexpected x: $x")
end

# single operator
for op in (:(𝒥𝓏), :(𝒥₊), :(𝒥₋), :(ℐ𝓏), :(ℐ₊), :(ℐ₋))
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
    return Operator(Matrix(c'), basis, basis)
    # return Operator(float.(c'), basis, basis)
end

𝒥₊ℐ₋(basis::AbstractVector{<:AtomState})  = ops(basis, (𝒥₊, ℐ₋))
𝒥₋ℐ₊(basis::AbstractVector{<:AtomState})  = ops(basis, (𝒥₋, ℐ₊))
𝒥₊²ℐ₋²(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥₊, 𝒥₊, ℐ₋, ℐ₋))
𝒥₋²ℐ₊²(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥₋, 𝒥₋, ℐ₊, ℐ₊))