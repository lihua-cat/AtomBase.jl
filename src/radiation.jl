@doc raw"""

    transitionME(J, I, F, MF, J', I', F', MF', k = 1)

return coefficient $c$ in $<JIFMF|P^k_q|J'I'F'MF'> = c \cdot <J||P^k||J'>$
"""
function transitionME(J1, I1, F1, MF1, J2, I2, F2, MF2, k::Int = 1)
    #   <FMF|T^k_q|F'MF'> = c1 * <F||T^k||F'>
    q = MF1 - MF2
    c1 = wigner_eckart(F1, F2, MF1, MF2, k, q)
    #   <JIF||T^k||J'I'F'> = c2 * <J||T^k||J>
    c2 = uncoup_T1(J1, I1, F1, J2, I2, F2, 1)
    return c1 * c2
end

@doc raw"""

    reducedME_M1(L, S, J, L', S', J')

no dimension value(/ħ) of $<LSJ||J + (gs - 1)S||L'S'J'>$.
"""
function reducedME_M1(L1, S1, J1, L2, S2, J2)
    gs = 2.00232
    #   <LSJ||(J+(gs-1)S)||L'S'J'> = <LSJ||J||L'S'J'> + (gs - 1)<LSJ||S||L'S'J'>
    rme_J = reduceME(J1, J2)
    #   <LSJ||S||L'S'J'> = c * <S||S||S>
    c = uncoup_T2(L1, S1, J1, L2, S2, J2, 1)
    rme_S = reduceME(S1, S2)
    return rme_J + (gs - 1) * c * rme_S
end

@doc raw"""

    reducedME_E1(L, S, J, L', S', J')

no dimension value(/a_0) of $<LSJ||r||L'S'J'>$.
"""
function reducedME_E1(L1, S1, J1, L2, S2, J2)
    #   <LSJ||r||L'S'J'> = c * <L||r||L>
    #   <L||r||L'> depends on ψ(r), radical wave function.
    c1 = uncoup_T1(L1, S1, J1, L2, S2, J2, 1)
    c2 = NaN    # need more information of ψ(r)
    c = c1 * c2
    return c
end

function relative_transition_intensity(
        b::Bra{T1, HyperfineStructureState{L1,S1,J1,I1}}, 
        k::Ket{T2, HyperfineStructureState{L2,S2,J2,I2}},
        order::String
        ) where {T1,T2,L1,S1,J1,I1,L2,S2,J2,I2}
    F1, F2 = b.s.F, k.s.F
    MF1, MF2 = b.s.MF, k.s.MF
    q = MF1 - MF2
    kk = parse(Int64, order[end])
    if abs(q) > kk
        c = 0
    else
        c1 = b.c * k.c * transitionME(J1, I1, F1, MF1, J2, I2, F2, MF2, kk)
        if order == "E1"
            c = c1 * reducedME_E1(L1, S1, J1, L2, S2, J2)
        elseif order == "M1"
            c = c1 * reducedME_M1(L1, S1, J1, L2, S2, J2)
        end
    end
    return c
end

function relative_transition_intensity(bv::BraVec, kv::KetVec, order::String)
    c = 0
    for b in bv, k in kv
        b.c * k.c ≈ 0 && continue
        c += relative_transition_intensity(b, k, order)
    end
    return c
end