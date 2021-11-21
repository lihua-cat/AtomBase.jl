@doc raw"""

    transitionME(J, I, F, MF, q)

compute transition matrix element $\langle \alpha JIFM_F|\mathrm{(J+S)}_q|\alpha' J' I' F' M_F'\rangle$.
"""
function transitionME(
        J::NTuple{2,T},
        I::NTuple{2,T}, 
        F::NTuple{2,T}, 
        MF::NTuple{2,T},
        q::T,
        k::Int = 1) where {T<:HalfInteger}
    #   <FMF|T^k_q|F'MF'> = c1 * <F||T^k||F'>
    c1 = wigner_eckart(F[1], F[2], MF[1], MF[2], q, k)
    #   <JIF||T^k||J'I'F'> = c2 * <J||T^k||J>
    c2 = uncoup_T1(J[1], I[1], F[1], J[2], I[2], F[2], 1)
    return c1 * c2
end

@doc raw"""

    reducedME_M1(L, S, J)

MD transition matrix element $\langle \alpha LSJ|\mathrm{(J+S)}_q|\alpha' L'S'J'\rangle$.
"""
function reducedME_M1(L, S, J)
    #   <LSJ||(J+S)||L'S'J'> = <LSJ||J||L'S'J'> + <LSJ||S||L'S'J'>
    rme_J = reduceME(J[1], J[2])
    #   <LSJ||S||L'S'J'> = c * <S||S||S>
    c = uncoup_T2(L[1], S[1], J[1], L[2], S[2], J[2], 1)
    rme_S = reduceME(S[1], S[2])
    if rme_J == 0
        return c * rme_S
    else
        return rme_J + c * rme_S
    end
end

@doc raw"""

    reducedME_E1(L, S, J)

(realtive) ED transition matrix element $\langle \alpha LSJ|\mathrm{r}_q|\alpha' L'S'J'\rangle$.
"""
function reducedME_E1(L, S, J)
    #   <LSJ||r||L'S'J'> = c * <L||r||L>
    #   <L||r||L'> depends on ψ(r), radical wave function.
    c = uncoup_T1k(L[1], S[1], J[1], L[2], S[2], J[2], 1)
    return c
end

function relative_transition_intensity(
        b::Bra{T1, HyperfineStructureState{L1,S1,J1,I1}}, 
        k::Ket{T2, HyperfineStructureState{L2,S2,J2,I2}}, 
        q::Int) where {T1,T2,L1,S1,J1,I1,L2,S2,J2,I2}
    L = HalfInt.(L1, L2)
    S = HalfInt.(S1, S2)
    J = HalfInt.(J1, J2)
    I = HalfInt.(I1, I2)
    F = (b.F, k.F)
    MF = (b.MF, k.MF)
    if (-1)^L1 == (-1)^L2
        c = transitionME(J, I, F, MF, q, 1) * reducedME_M1(L, S, J)
    else
        c = transitionME(J, I, F, MF, q, 1) * reducedME_E1(L, S, J)
    end
    return c
end

function relative_transition_intensity(bv::BraVec, kv::KetVec, q::Int)
    c = 0
    for b in bv, k in kv
        c += relative_transition_intensity(b, k, q)
    end
    return c
end