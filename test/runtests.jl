using Test
using AtomBase
using LinearAlgebra

signshow(m) =
    if m > 0
        '+' * string(Int(m))
    elseif m == 0
        ' ' * string(Int(m))
    else
        string(Int(m))
    end

L = 2
S = 3 // 2
In = 3 // 2
for J = (L+S):-1:abs(L-S)
    hfs_c = basis_hfs(L, S, J, In; couple = true)
    hfs_uc = basis_hfs(L, S, J, In; couple = false)
    @testset "J=$J I=$In Basis          " begin
        @test hfs_c.F == (J+In):-1:abs(J - In)
        @test hfs_uc.MF == (J+In):-1:(-(J + In))
        @test length.(hfs_c.B) == 2hfs_c.F .+ 1
        @test length.(hfs_uc.B1) == length.(hfs_uc.B2)
        @test sum(length.(hfs_c.B)) == sum(length.(hfs_uc.B1))
    end
    for MF = (J+In):-1:(-(J+In))
        basis = basis_get(hfs_uc, :B1, :MF, MF)
        basis2 = basis_get(hfs_uc, :B2, :MF, MF)
        c = normalize(rand(length(basis)))
        k = KetVec(c, basis)
        b = k'
        @testset "J=$J I=$In MF=$(signshow(MF)) Ket Bra  " begin
            @test k' * k ≈ 1
            @test b * k ≈ 1
        end
        @testset "J=$J I=$In MF=$(signshow(MF)) Operator " begin
            Jz = 𝒥𝓏(basis)
            Iz = ℐ𝓏(basis)
            @test Jz.c == Diagonal([s.MJ for s in basis])
            @test Iz.c == Diagonal([s.MI for s in basis])
            val, vec = diagonal(Jz)
            @test val == [s.MJ for s in basis]
            J₊I₋ = 𝒥₊ℐ₋(basis)
            J₋I₊ = 𝒥₋ℐ₊(basis)
            if J >= In
                val1, vec1 = eigen(J₊I₋.c + J₋I₊.c)
            else
                val1, vec1 = eigen(J₊I₋.c + J₋I₊.c; sortby = x -> -x)
            end
            val2, vec2 = diagonal(J₊I₋ + J₋I₊)
            @test val1 == val2
        end
        @testset "J=$J I=$In MF=$(signshow(MF)) Transofrm" begin
            kt = basistransform(k, basis2)
            ktt = basistransform(kt, basis)
            @test k.c ≈ ktt.c && k.s == ktt.s
        end
    end
end