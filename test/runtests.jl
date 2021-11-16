using Test
using AtomBase
using LinearAlgebra

signshow(m) = if m > 0
    '+' * string(Int(m))
elseif m == 0
    ' ' * string(Int(m))
else
    string(Int(m))
end

L = 2
S = 3//2
In = 3//2
for J in (L + S):-1:abs(L - S)
    hfs_c = basis_hfs(L, S, J, In; couple=true)
    hfs_uc = basis_hfs(L, S, J, In; couple=false)
    @testset "J=$J I=$In Basis          " begin
        @test hfs_c.F == (J + In):-1:abs(J - In)
        @test hfs_uc.MF == (J + In):-1:(-(J + In))
        @test length.(hfs_c.B) == 2hfs_c.F .+ 1
        @test length.(hfs_uc.B1) == length.(hfs_uc.B2)
        @test sum(length.(hfs_c.B)) == sum(length.(hfs_uc.B1))
        for MF in (J + In):-1:(-(J + In))
            basis = basis_get(hfs_uc, :B1, :MF, MF)
            @test basis.MJ + basis.MI == fill(MF, length(basis))
        end
    end
    for MF in (J + In):-1:(-(J + In))
        basis = basis_get(hfs_uc, :B1, :MF, MF)
        basis2 = basis_get(hfs_uc, :B2, :MF, MF)
        c = normalize(rand(length(basis)))
        k = Ket(c, basis)
        b = Bra(c', basis')
        @testset "J=$J I=$In MF=$(signshow(MF)) Ket Bra  " begin
            @test k' * k ≈ 1
            @test k' == b
            @test b * k ≈ 1
        end
        @testset "J=$J I=$In MF=$(signshow(MF)) Operator " begin
            Jz = 𝒥𝓏(basis)
            Iz = ℐ𝓏(basis)
            @test Jz.c == Diagonal([M for M in basis.MJ])
            @test Iz.c == Diagonal([M for M in basis.MI])
            val, vec = eigen(Jz)
            @test val == [M for M in basis.MJ]
            J₊I₋ = 𝒥₊ℐ₋(basis)
            J₋I₊ = 𝒥₋ℐ₊(basis)
            if J >= In
                val1, vec1 = eigen(J₊I₋.c + J₋I₊.c)
            else
                val1, vec1 = eigen(J₊I₋.c + J₋I₊.c; sortby=x -> -x)
            end
            val2, vec2 = eigen(J₊I₋ + J₋I₊)
            @test val1 == val2
            @test vec1 == vec2.c
        end
        @testset "J=$J I=$In MF=$(signshow(MF)) Transofrm" begin
            kt = transform(k, basis2)
            ktt = transform(kt, basis)
            @test k.c ≈ ktt.c && k.s == ktt.s
        end
    end
end