abstract type Dirac{T,S} end

struct Ket{T<:Number,S<:AtomState} <: Dirac{T,S}
    c::T
    s::S
end
Ket(s::AtomState) = Ket(1, s)

struct Bra{T<:Number,S<:AtomState} <: Dirac{T,S}
    c::T
    s::S
end
Bra(s::AtomState) = Bra(1, s)

struct KetVec{T<:Number,S<:AtomState} <: AbstractVector{Ket{T,S}}
    c::Vector{T}
    s::Vector{S}
end

const BraVec{T,S} = Adjoint{Bra{T,S},KetVec{T,S}}

struct Op{T<:Number,S<:AtomState}
    c::T
    ks::S
    bs::S
end

struct Operator{T<:Number,S<:AtomState} <: AbstractMatrix{Op{T,S}}
    c::Matrix{T}
    ks::Vector{S}
    bs::Vector{S}
end

## Index
Base.size(k::KetVec) = size(k.c)
Base.getindex(k::KetVec, i) = Ket(k.c[i], k.s[i])
Base.size(op::Operator) = size(op.c)
Base.getindex(op::Operator, i, j) = Op(op.c[i, j], op.ks[i], op.bs[j])

## LinearAlgebra
#  adjoint
adjoint(k::Ket) = Bra(k.c', k.s)
adjoint(b::Bra) = Ket(b.c', b.s)
adjoint(op::Op) = Op(op.c', op.bs, op.ks)
# zero(d)
zero(d::D) where {D<:Dirac} = D(zero(d.c), d.s)
zero(op::Op) = Op(zero(op.c), op.ks, op.bs)
#  ket + ket, bra + bra
+(d1::D, d2::D) where {D<:Dirac} =
    d1.s == d2.s ? D(d1.c + d2.c, d1.s) : error("ket or bra +")
+(k1::KetVec, k2::KetVec) = k1.s == k2.s ? KetVec(k1.c + k2.c, k1.s) : error("KetVec +")
+(b1::BraVec, b2::BraVec) = (b1' + b2')'
#  op + op
+(op1::Op{T1,S}, op2::Op{T2,S}) where {T1,T2,S} =
    op1.ks == op2.ks && op1.bs == op2.bs ? Op(op1.c + op2.c, op1.ks, op2.ks) : error("op +")
+(op1::Operator{T1,S}, op2::Operator{T2,S}) where {T1,T2,S} =
    op1.ks == op2.ks && op1.bs == op2.bs ? Operator(op1.c + op2.c, op1.ks, op2.ks) : error("op +")
# -ket or bra
-(d::D) where {D<:Dirac} = D(-d.c, d.s)
-(d1::D, d2::D) where {D<:Dirac} = d1 + (-d2)
-(k::KetVec) = KetVec(-k.c, k.s)
-(k1::KetVec, k2::KetVec) = k1 + (-k2)
# -op
-(op::Op) = Op(-op.c, op.ks, op.bs)
-(op1::Op, op2::Op) = op1 + (-op2)
-(op::Operator) = Operator(-op.c, op.ks, op.bs)
-(op1::Operator, op2::Operator) = op1 + (-op2)
#  x * dirac, dirac * x
*(x::Number, d::D) where {D<:Dirac} = x == 0 ? x : D(x * d.c, d.s)
*(d::D, x::Number) where {D<:Dirac} = *(x, d)
#  bra * ket
*(b::Bra{T1,S}, k::Ket{T2,S}) where {T1,T2,S} = b.s == k.s ? b.c * k.c : zero(b.c * k.c)
#  ket * bra = op
*(k::Ket{T1,S}, b::Bra{T2,S}) where {T1,T2,S} =
    Op{promote_type(T1, T2),S}(k.c * b.c, k.s, b.s)
#  op * ket
*(op::Op{T1,S}, k::Ket{T2,S}) where {T1,T2,S} = Ket(Bra(op.c, op.bs) * k, op.ks)
#  x * op, op * x
*(x::Number, op::Op) = Op(x * op.c, op.ks, op.bs)
*(op::Op, x::Number) = *(x, op)
*(x::Number, op::Operator) = Operator(x * op.c, op.ks, op.bs)
*(op::Operator, x::Number) = *(x, op)
## Quantum Mechanics
