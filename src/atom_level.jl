abstract type EnergyLevel end

struct FineLevel <: EnergyLevel
    state::FineState
    hfc::HyperfineConstant
    E::Wavenumber   # energy level
end