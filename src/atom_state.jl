"""
quantum number notation
- `L` electron orbital quantum number
- `S` electron spin quantum number
- `J` L + S, total quantum number of fine structure
- `I` nuclear spin quantum number
- `F` J + I, total quantum number of hyperfine sturcture 
"""
abstract type State end
abstract type SubState <: State end

## Fine Structure
struct FineState <: State
    L::Int
    S::HalfInt
    J::HalfInt
end

struct CoupledSubFineState <: SubState
    L::Int
    S::HalfInt
    J::HalfInt
    MJ::HalfInt
end

struct UncoupledSubFineState <: SubState
    L::Int
    ML::Int
    S::HalfInt
    MS::HalfInt
end

## Hyperfine structure
struct HyperfineState <: State
    L::Int
    S::HalfInt
    J::HalfInt
    I::HalfInt
    F::HalfInt
end

struct CoupledSubHyperfineState <: SubState
    L::Int
    S::HalfInt
    J::HalfInt
    I::HalfInt
    F::HalfInt
    MF::HalfInt
end

struct UncoupledSubHyperfineState <: SubState
    L::Int
    S::HalfInt
    J::HalfInt
    MJ::HalfInt
    I::HalfInt
    MI::HalfInt
end