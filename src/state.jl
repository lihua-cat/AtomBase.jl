abstract type AtomState end

struct HyperfineStructureState{L,S,J,I} <: AtomState
    F::HalfInt
    MF::HalfInt
    function HyperfineStructureState{L,S,J,I}(F, MF) where {L,S,J,I}
        @check_args HyperfineStructureState F >= 0
        @check_args HyperfineStructureState MF in F:-1:(-F)
        new{L,S,J,I}(F, MF)
    end
end

struct UncoupledHyperfineStructureState{L,S,J,I} <: AtomState
    MJ::HalfInt
    MI::HalfInt
    function UncoupledHyperfineStructureState{L,S,J,I}(MJ, MI) where {L,S,J,I}
        @check_args UncoupledHyperfineStructureState J >= 0
        @check_args UncoupledHyperfineStructureState MJ in J:-1:(-J)
        @check_args UncoupledHyperfineStructureState I > 0
        @check_args UncoupledHyperfineStructureState MI in I:-1:(-I)
        new{L,S,J,I}(MJ, MI)
    end
end