"""
quantum number notation
- `L` electron orbital quantum number
- `S` electron spin quantum number
- `J` L + S, total quantum number of fine structure
- `I` nuclear spin quantum number
- `F` J + I, total quantum number of hyperfine sturcture 
"""
abstract type AtomState end

struct FineStructure <: AtomState
    L::Int
    S::HalfInt
    J::HalfInt
    function FineStructure(L, S, J)
        if all((L, S, J) .>= 0)
            if J in abs(L-S):abs(L+S)
                new(L, S, J)
            else
                error("L, S, J not match: $L, $S, $J")
            end
        else
            error("L: $L, S: $S, J: $J must be non-negative!!")
        end
    end
end

struct HyperfineStructure <: AtomState
    L::Int
    S::HalfInt
    J::HalfInt
    I::HalfInt
    F::HalfInt
    function FineStructure(L, S, J, I, F)
        if all((L, S, J, I, F) .>= 0)
            if J in abs(L-S):abs(L+S)
                if F in abs(J-I):abs(J+I)
                    new(L, S, J)
                else
                    error("J, I, F not match: $J, $I, $F")
                end
            else
                error("L, S, J not match: $L, $S, $J")
            end
        else
            error("L: $L, S: $S, J: $J, I: $I, F: $F must be non-negative!!")
        end
    end
end