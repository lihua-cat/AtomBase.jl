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


function hamiltonian_hfs(basis::AbstractVector{UncoupledHyperfineStructureState{L,S,J,I}}, A, B) where {L,S,J,I}
    Jz = 𝐉𝐳(basis)
    Iz = 𝐈𝐳(basis)
    J₊I₋ = 𝐉₊𝐈₋(basis)
    J₋I₊ = 𝐉₋𝐈₊(basis)
    JzJ₊IzI₋ = 𝐉𝐳𝐉₊𝐈𝐳𝐈₋(basis)
    JzJ₊I₋Iz = 𝐉𝐳𝐉₊𝐈₋𝐈𝐳(basis)
    J₊JzIzI₋ = 𝐉₊𝐉𝐳𝐈𝐳𝐈₋(basis)
    J₊JzI₋Iz = 𝐉₊𝐉𝐳𝐈₋𝐈𝐳(basis)
    JzJ₋IzI₊ = 𝐉𝐳𝐉₋𝐈𝐳𝐈₊(basis)
    JzJ₋I₊Iz = 𝐉𝐳𝐉₋𝐈₊𝐈𝐳(basis)
    J₋JzIzI₊ = 𝐉₋𝐉𝐳𝐈𝐳𝐈₊(basis)
    J₋JzI₊Iz = 𝐉₋𝐉𝐳𝐈₊𝐈𝐳(basis)
    J₊²I₋² = 𝐉₊²𝐈₋²(basis)
    J₋²I₊² = 𝐉₋²𝐈₊²(basis)
    E = Operator(Matrix{Int}(𝐼(length(basis))), basis, basis)
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