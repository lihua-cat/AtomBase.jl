function basis_hfs(L::Int, S::HalfInt, J::HalfInt, I::HalfInt; couple::Bool)
    if couple
        F = collect((J+I):-1:abs(J - I))
        num = length(F)
        B = Vector{Vector{HyperfineStructureState{L,S,J,I}}}(undef, num)
        for i = 1:num
            B[i] = [HyperfineStructureState{L,S,J,I}(F[i], MF) for MF in collect(F[i]:-1:(-F[i]))]
        end
        data = ndsparse((F = F,), (basis = B,))
    else
        MF = collect((J+I):-1:(-(J + I)))
        num = length(MF)
        B1 = Vector{Vector{UncoupledHyperfineStructureState{L,S,J,I}}}(undef, num)
        B2 = Vector{Vector{HyperfineStructureState{L,S,J,I}}}(undef, num)
        for i = 1:num
            if J >= I
                MJ = collect(max(-J, MF[i] - I):min(J, MF[i] + I))
                MI = MF[i] .- MJ
            else
                MI = collect(max(-I, MF[i] - J):min(I, MF[i] + J))
                MJ = MF[i] .- MI
            end
            B1[i] = [UncoupledHyperfineStructureState{L,S,J,I}(MJ[j], MI[j]) for j = 1:length(MJ)]
            B2[i] = [HyperfineStructureState{L,S,J,I}(F, MF[i]) for F in collect((J+I):-1:max(abs(J - I), abs(MF[i])))]
        end
        data = ndsparse((MF = MF,), (basis1 = B1, basis2 = B2))
    end
    return data
end

function basis_hfs(L, S, J, I; couple::Bool)
    return basis_hfs(Int(L), HalfInt(S), HalfInt(J), HalfInt(I); couple = couple)
end