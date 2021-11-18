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
    c = eltype(c) == Real ? float.(c) : c
    return Operator(Matrix(c'), basis, basis)
    # return Operator(float.(c'), basis, basis)
end

𝒥₊ℐ₋(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥₊, ℐ₋))
𝒥₋ℐ₊(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥₋, ℐ₊))

𝒥𝓏𝒥₊ℐ𝓏ℐ₋(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥𝓏, 𝒥₊, ℐ𝓏, ℐ₋))
𝒥𝓏𝒥₊ℐ₋ℐ𝓏(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥𝓏, 𝒥₊, ℐ₋, ℐ𝓏))
𝒥₊𝒥𝓏ℐ𝓏ℐ₋(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥₊, 𝒥𝓏, ℐ𝓏, ℐ₋))
𝒥₊𝒥𝓏ℐ₋ℐ𝓏(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥₊, 𝒥𝓏, ℐ₋, ℐ𝓏))

𝒥𝓏𝒥₋ℐ𝓏ℐ₊(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥𝓏, 𝒥₋, ℐ𝓏, ℐ₊))
𝒥𝓏𝒥₋ℐ₊ℐ𝓏(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥𝓏, 𝒥₋, ℐ₊, ℐ𝓏))
𝒥₋𝒥𝓏ℐ𝓏ℐ₊(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥₋, 𝒥𝓏, ℐ𝓏, ℐ₊))
𝒥₋𝒥𝓏ℐ₊ℐ𝓏(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥₋, 𝒥𝓏, ℐ₊, ℐ𝓏))

𝒥₊²ℐ₋²(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥₊, 𝒥₊, ℐ₋, ℐ₋))
𝒥₋²ℐ₊²(basis::AbstractVector{<:AtomState}) = ops(basis, (𝒥₋, 𝒥₋, ℐ₊, ℐ₊))


function hamiltonian_hfs(basis::AbstractVector{UncoupledHyperfineStructureState{L,S,J,I}}, A, B) where {L,S,J,I}
    Jz = 𝒥𝓏(basis)
    Iz = ℐ𝓏(basis)
    J₊I₋ = 𝒥₊ℐ₋(basis)
    J₋I₊ = 𝒥₋ℐ₊(basis)
    JzJ₊IzI₋ = 𝒥𝓏𝒥₊ℐ𝓏ℐ₋(basis)
    JzJ₊I₋Iz = 𝒥𝓏𝒥₊ℐ₋ℐ𝓏(basis)
    J₊JzIzI₋ = 𝒥₊𝒥𝓏ℐ𝓏ℐ₋(basis)
    J₊JzI₋Iz = 𝒥₊𝒥𝓏ℐ₋ℐ𝓏(basis)
    JzJ₋IzI₊ = 𝒥𝓏𝒥₋ℐ𝓏ℐ₊(basis)
    JzJ₋I₊Iz = 𝒥𝓏𝒥₋ℐ₊ℐ𝓏(basis)
    J₋JzIzI₊ = 𝒥₋𝒥𝓏ℐ𝓏ℐ₊(basis)
    J₋JzI₊Iz = 𝒥₋𝒥𝓏ℐ₊ℐ𝓏(basis)
    J₊²I₋² = 𝒥₊²ℐ₋²(basis)
    J₋²I₊² = 𝒥₋²ℐ₊²(basis)
    E = Operator(Matrix{Int}(𝐼, length(basis), length(basis)), basis, basis)
    h_md = A * (Jz * Iz + (J₊I₋ + J₋I₊) / 2)
    if J == 1 / 2 && I >= 1
        h_hfs = h_md
    else
        c1 = B / J / (2J - 1) / 2I / (2I - 1)
        c2 = 3 * Jz * Jz - J * (J + 1) * E
        c3 = 3 * Iz * Iz - I * (I + 1) * E
        c4 = JzJ₊IzI₋ + JzJ₊I₋Iz + J₊JzIzI₋ + J₊JzI₋Iz +
             JzJ₋IzI₊ + JzJ₋I₊Iz + J₋JzIzI₊ + J₋JzI₊Iz + J₋²I₊² + J₊²I₋²
        h_eq = c1 * (1 // 2 * c2 * c3 + 3 // 4 * c4)
        h_hfs = h_md + h_eq
    end
    return h_hfs
end