abstract type EnergyLevel end

struct FineLevel <: EnergyLevel
    state::FineStructure
    hfc::HyperfineConstant
    E::Wavenumber   # energy level
end