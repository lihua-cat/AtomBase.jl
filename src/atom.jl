"""
Atom
- `M::Mass`                         relative atomic mass, default in unit `u`
- `gI::Float64`                     nuclear Lande g-factor
- `I::HalfInt`                      nuclear spin quantum number, either an integer or a half-integer
- `ground::FineLevel`               only one ground state (fine stucture)
- `excited::Vector{FineLevel}`      several excited states (fine stucture)
"""
struct Atom
    M::Mass                         # relative atomic mass, default in unit `u`
    gI::Float64                     # nuclear Lande g-factor
    I::HalfInt                      # nuclear spin quantum number, either an integer or a half-integer
    ground::FineLevel               # only one ground state (fine stucture)
    excited::Vector{FineLevel}      # several excited states (fine stucture)
end