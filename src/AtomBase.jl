module AtomBase

using DataFrames
using Printf
using HalfIntegers, RationalRoots
using WignerSymbols
using LinearAlgebra
import LinearAlgebra: adjoint, Adjoint, I as 𝐼
import Base: +, -, *, /, zero


include("utils.jl")

export AtomState
export HyperfineStructureState, UncoupledHyperfineStructureState
include("state.jl")

export basis_hfs, basis_get
include("basis.jl")

export Dirac, Ket, Bra, KetVec, BraVec, Op, Operator
include("dirac.jl")

export 𝒥𝓏, 𝒥₊, 𝒥₋
export ℐ𝓏, ℐ₊, ℐ₋
export 𝒥₊ℐ₋, 𝒥₋ℐ₊
export 𝒥₊²ℐ₋², 𝒥₋²ℐ₊²
export hamiltonian_hfs
include("operator.jl")

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
