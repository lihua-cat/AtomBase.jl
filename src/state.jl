abstract type AtomState end

struct HyperfineStructureState{L,S,J,I} <: AtomState
    F::HalfInt
    MF::HalfInt
    function HyperfineStructureState{L,S,J,I}(F, MF) where {L,S,J,I}
        if F >= 0 && MF in F:-1:(-F)
            new{L,S,J,I}(F, MF)
        else
            error("Invalid angular momentum quantum number")
        end
    end
end

struct UncoupledHyperfineStructureState{L,S,J,I} <: AtomState
    MJ::HalfInt
    MI::HalfInt
    function UncoupledHyperfineStructureState{L,S,J,I}(MJ, MI) where {L,S,J,I}
        if J > 0 && MJ in J:-1:(-J) && I > 0 && MI in I:-1:(-I)
            new{L,S,J,I}(MJ, MI)
        else
            error("Invalid angular momentum quantum number")
        end
    end
end