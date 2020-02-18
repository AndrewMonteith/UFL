abstract type Foo end 

struct Concrete1 <: Foo 
    a::Int

    function Concrete1(a::Int)
        println("CTOR 1 called")
        new(a)
    end
end 

struct Concrete2 <: Foo 
    a::Int 

    function Concrete2(b::Int)
        println("CTOR 2 called")
        new(b)
    end
end


function do_something()

    xs::Vector{Foo} = [Concrete1(1), Concrete2(3)]
    new_vec::Vector{Foo} = []

    for x ∈ xs
        y = typeof(x)(x.a)

        push!(new_vec, y)
    end

    new_vec
end


do_something()