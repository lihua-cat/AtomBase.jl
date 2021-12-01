module AtomBase

import PhysicalConstants.CODATA2018: h as ℎ, μ_B as 𝜇B, μ_0 as 𝜇0, ε_0 as 𝜀₀, e as 𝑒, a_0 as 𝑎₀

using DataFrames
using Printf
using HalfIntegers, RationalRoots
using WignerSymbols
using LinearAlgebra
import LinearAlgebra: adjoint, Adjoint, I as 𝐼
import Base: +, -, *, /, zero

using UsefulFunctions

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

export wigner_eckart, reduceME, uncoup_T1, uncoup_T2
include("tensor.jl")

export transitionME, reducedME_M1, reducedME_E1, relative_transition_intensity
include("radiation.jl")

include("show.jl")


end
