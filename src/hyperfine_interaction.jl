"""
HyperfineConstant
- `A::Wavenumber{<:Real}`   magnetic dipole constant
- `B::Wavenumber{<:Real}`   electric quadrupole constant
"""
struct HyperfineConstant
    A::Wavenumber{<:Real}   # magnetic dipole constant
    B::Wavenumber{<:Real}   # electric quadrupole constant
end