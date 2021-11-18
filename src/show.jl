function string_coe(c)
    if c isa AbstractFloat
        if abs(c) < 1e-6
            c_str = @sprintf "%9.3e" c
        else
            c_str = @sprintf "%9.6f" c
        end
    elseif c isa Integer
        c_str = @sprintf "%-2d" c
    else
        c_str = repr(c, context = :compact => true)
    end
end

## AtomState
function Base.show(io::IO, s::AtomState)
    for p in propertynames(s)
        v = getproperty(s, p)
        if p in (:MJ, :MI, :MF, :MS)
            v > 0 && print(io, "+$v")
            v < 0 && print(io, "$v")
            v == 0 && print(io, " $v")
        else
            print(io, "$v")
        end
        p == propertynames(s)[end] || print(io, ",")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", s::AtomState)
    println(io, summary(s))
    for p in propertynames(s)
        v = getproperty(s, p)
        if p in (:MJ, :MI, :MF, :MS)
            v > 0 && print(io, "$p=+$v")
            v < 0 && print(io, "$p=$v")
            v == 0 && print(io, "$p= $v")
        else
            print(io, "$p=$v")
        end
        p == propertynames(s)[end] || print(io, ",")
    end
    return nothing
end

## Dirac
function Base.show(io::IO, d::Dirac)
    (;c, s) = d
    c_str = string_coe(c)
    printstyled(io, c_str, color = :blue)
    d isa Ket ? print(io, "|") : print(io, "<")
    print(io, s)
    d isa Ket ? print(io, ">") : print(io, "|")
    return nothing
end

function Base.show(io::IO, op::Op)
    (;c, ks, bs) = op
    c_str = string_coe(c)
    printstyled(io, c_str, color = :blue)
    print(io, "|")
    print(io, ks)
    print(io, ">")
    print(io, "<")
    print(io, bs)
    print(io, "|")
    return nothing
end

function Base.show(io::IO, mime::MIME"text/plain", d::Dirac)
    println(io, summary(d))
    print(io, ' ')
    print(io, d)
    return nothing
end

function Base.show(io::IO, mime::MIME"text/plain", op::Op)
    println(io, summary(op))
    print(io, ' ')
    print(io, op)
    return nothing
end
