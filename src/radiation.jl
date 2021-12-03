@doc raw"""

    relative_transitionME(J, I, F, MF, J', I', F', MF', k)

return coefficient $C$ in 
$\bra{\alpha JIFM_F}\mathbf{P}^k_q\ket{\alpha' J'I'F'M'_F} = 
C \cdot \bra{\alpha J}|\mathbf{P}^k|\ket{\alpha' J'}$
"""
function relative_transitionME(J1, I1, F1, MF1, J2, I2, F2, MF2, k)
    #   <FMF|T^k_q|F'MF'> = c1 * <F||T^k||F'>
    q = MF1 - MF2
    c1 = wigner_eckart(F1, F2, MF1, MF2, k, q)
    #   <JIF||T^k||J'I'F'> = c2 * <J||T^k||J>
    c2 = uncoup_T1(J1, I1, F1, J2, I2, F2, k)
    return c1 * c2
end

@doc raw"""

    reducedME_E1(L, S, J, L', S', J')

value of $\bra{LSJ}|\mathbf{P}^{(1)}|\ket{L'S'J'}$.
```math
\mathbf{P}^{(1)} = -\frac{e}{\sqrt{4\pi\epsilon_0}}\mathbf{r}^{(1)}
```
"""
function reducedME_E1(L1, S1, J1, L2, S2, J2)
    #   <LSJ||r||L'S'J'> = c * <L||r||L>
    #   <L||r||L'> depends on ψ(r), radical wave function.
    c = uncoup_T1(L1, S1, J1, L2, S2, J2, 1)
    rme = NaN # <L||r||L'> / a_0(Bohr radius)
    ME = c * rme
    C = -𝑒 * 𝑎₀ / √(4π * 𝜀₀)    # from cgs unit to SI unit
    ME *= C
    return ME
end

@doc raw"""

    reducedME_M1(L, S, J, L', S', J')

value of $\bra{LSJ}|\mathbf{P}^{(1)}|\ket{L'S'J'}$.
```math
\mathbf{P}^{(1)} = -\sqrt{\frac{\mu_0}{4\pi}}\mu_B\left[\mathbf{J}^{(1)}+(g_s-1)\mathbf{S}^{(1)}\right]
```
"""
function reducedME_M1(L1, S1, J1, L2, S2, J2, gs = gS)
    #   <LSJ||(J+(gs-1)S)||L'S'J'> = <LSJ||J||L'S'J'> + (gs - 1)<LSJ||S||L'S'J'>
    rme_J = reducedME(J1, J2)
    #   <LSJ||S||L'S'J'> = c * <S||S||S>
    c = uncoup_T2(L1, S1, J1, L2, S2, J2, 1)
    rme_S = reducedME(S1, S2)
    ME = rme_J + (gs - 1) * c * rme_S   # / ħ
    C = -√(𝜇0 / 4π) * 𝜇B    # from cgs unit to SI unit
    ME *= C
    return ME
end

function relative_transitionME(
        b::Bra{T1, HyperfineStructureState{L1,S1,J1,I1}}, 
        k::Ket{T2, HyperfineStructureState{L2,S2,J2,I2}},
        order::Int
        ) where {T1,T2,L1,S1,J1,I1,L2,S2,J2,I2}
    F1, F2 = b.s.F, k.s.F
    MF1, MF2 = b.s.MF, k.s.MF
    q = MF1 - MF2
    if abs(q) > order
        c = 0
    else
        c = b.c * k.c * relative_transitionME(J1, I1, F1, MF1, J2, I2, F2, MF2, order)
    end
    return c
end

@doc raw"""

    relative_transitionME(ψⱼ::BraVec, ψᵢ::KetVec, k::Int)

return coefficient $c$ in 
$\bra{\psi_j}\mathbf{P}^k\ket{\psi_i} = c \cdot \bra{\alpha J}|\mathbf{P}^k|\ket{\alpha' J'}$. 
$\ket{\psi}$ is a ket in $\ket{\alpha J I FM_F}$ basis with the same quantum number $M_F$.
```math
\ket{\psi_j} = \sum_F c^j_F \ket{\alpha J I FM_F},\ \ket{\psi_i} = \sum_{F'} c^i_{F'} \ket{\alpha J' I' F'M'_F}
```
```math
c = \sum_{F, F'} {(c^j_F)}^\ast c^i_{F'} C(J, I, F, M_F, J', I', F', M'_F)
```
"""
function relative_transitionME(bv::BraVec, kv::KetVec, order::Int)
    c = 0
    for b in bv, k in kv
        b.c * k.c ≈ 0 && continue
        c += relative_transitionME(b, k, order)
    end
    return c
end

@doc raw"""

    aᵢⱼ(k, s)

Convert the transition matrix element between two state $j$ and $i$ to the Einstein spontaneous 
emission transition probability rate $a_{ij}$ (in unit $s^{-1}$)
"""
aᵢⱼ(k, s) = 64π^4 / 3ℎ * k^3 * s |> u"s^(-1)"

@doc raw"""

    einsteinA(k, L, S, J, I, F, L', S', J', I', F', term)

The Einstein A coefficient: spontaneous emission transition probability rate 
from a *specific* state $i$ of upper upper level to any state *j* of lower level.
```math
A = \frac{64\pi^4 k^3}{3h(2F'+1)} \bra{LSJIF}|\textrm{term of radiation}|\ket{L'S'J'I'F'}
```
"""
function einsteinA(k, L1, S1, J1, I1, F1, L2, S2, J2, I2, F2, term)
    order = parse(Int, term[end])
    c = uncoup_T1(J1, I1, F1, J2, I2, F2, order) |> Float64
    if term == "E1"
        rme = reducedME_E1(L1, S1, J1, L2, S2, J2)
    elseif term == "M1"
        rme = reducedME_M1(L1, S1, J1, L2, S2, J2)
    else
        error("term: $term is undefined")
    end
    g = 2F2 + 1
    s = (c * rme)^2
    gA = aᵢⱼ(k, s)
    A = gA / g
    return A
end

@doc raw"""

    σᵢⱼ(k, a, lineshape)

Convert the transition probability rate $a_{ij}$ to transition cross section 
$\sigma_{ij}$ (in unit $cm^2$) in wavenumber `k`.
```math
\sigma_{ij} = \frac{a_{ij}}{8\pi k^2} * f(k)
```
where $f(k)$ is line bradening value at wavenumber $k$. Gain between these two levels is given by
```math
g = \sum_{ij} \sigma_{ij}\left(N_j - N_i \right) = \sum_{ij} \sigma_{ij} \Delta N,\ \Delta N = \frac{N_u}{g_u} - \frac{N_l}{g_l}
```
"""
σᵢⱼ(k, a, lineshape) = 1 / k^2 / 8π * a * lineshape |> u"cm^2"
