
# Base.iterate(i::Int32) = (0, 0)
# function Base.iterate(i::Int32, state::Int32)
#     i === state && return nothing 

#     (state += 1, state)
# end

# function foo() 
#     for x in 10 
#         println(x)
#     end
# end

# foo()


# struct Foobar 
#     x::Int 
# end 

# Base.:(==)(x::Foobar, y::Foobar) = abs(x.x) === abs(y.x) #Base.isequal(x, y)

# f1, f2 = Foobar(-1), Foobar(1)
# f3, f4 = Foobar(2), Foobar(1)

# println(f1 == f2)
# println(f3 != f4)

struct Foo end 

function Base.:+(f::Foo, x)
    println("did some + stuff")
    f 
end

f = Foo()
println(2 + f)

# function Base.:-(f::Foo)
#     println("got unary minus")
#     f 
# end


# f = Foo()

# println(-f)

# struct Foobar 
#     function Foobar(a::Int)
#         if a === 0 
#             Foo()
#         else
#             new()
#         end
#     end
# end


# f = Foobar(0)

# println(f)
# println(typeof(f))