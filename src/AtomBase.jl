module AtomBase

import PhysicalConstants.CODATA2018: h as ℎ, μ_B as 𝜇B, μ_0 as 𝜇0, ε_0 as 𝜀₀, e as 𝑒, a_0 as 𝑎₀

using DataFrames
using Unitful
using Printf
using HalfIntegers, RationalRoots
using WignerSymbols
using LinearAlgebra
import LinearAlgebra: adjoint, Adjoint, I as 𝐼
import Base: +, -, *, /, zero

using UsefulFunctions

export gS
const gS = 1.00115965218085 * 2 #https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.030801

export AtomState
export HyperfineStructureState, UncoupledHyperfineStructureState
include("state.jl")

export basis_hfs
include("basis.jl")

export Dirac, Ket, Bra, KetVec, BraVec, Op, Operator
include("dirac.jl")

# export 𝐉𝐳, 𝐉₊, 𝐉₋
# export 𝐈𝐳, 𝐈₊, 𝐈₋
# export 𝐉₊𝐈₋, 𝐉₋𝐈₊
# export 𝐉₊²𝐈₋², 𝐉₋²𝐈₊²
include("operator.jl")

export lande
export hamiltonian_zeeman, hamiltonian_hfs, hamiltonian_total
include("hamiltonian.jl")

export diagonal
include("perturbation.jl")

export basistransform
include("basis transform.jl")

export wigner_eckart, reducedME, uncoup_T1, uncoup_T2
include("tensor.jl")

export relative_transitionME, reducedME_M1, reducedME_E1
export aᵢⱼ, einsteinA, σᵢⱼ
include("radiation.jl")

include("show.jl")


end
