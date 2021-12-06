function lande(l, s, j; gl = 1, gs = gS)
    j² = j * (j + 1) 
    l² = l * (l + 1)
    s² = s * (s + 1)
    gj = gl * (j² + l² - s²) / 2j² + gs * (j² - l² + s²) / 2j²
    return gj
end

function hamiltonian_zeeman(basis::AbstractVector{UncoupledHyperfineStructureState{L,S,J,I}},
    B; μB, μN, gl, gs, gI) where {L,S,J,I}
    gJ = lande(L, S, J, gl = gl, gs = gs)
    Jz = 𝐉𝐳(basis)
    Iz = 𝐈𝐳(basis)
    μ = -μB * gJ * Jz + μN * gI * Iz
    h = -μ * B
    return h
end

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
    if J == 1 / 2 || I < 1
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

function hamiltonian_total(basis, A, B, BF; μB, μN, gl, gs, gI)
    h_zeeman = hamiltonian_zeeman(basis, BF; μB = μB, μN = μN, gl = gl, gs = gs, gI = gI)
    h_hfs = hamiltonian_hfs(basis, A, B)
    h_total = h_zeeman + h_hfs
    return h_total
end