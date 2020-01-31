abstract type Foo end 
struct FooBar <: Foo end 

Base.getindex(f::Foo, x...) = x[2]

function Base.getindex(f::FooBar, x...)
    if x[1] === 0 
        return 100 
    else
        invoke(Base.getindex, Tuple{Foo, Vararg}, f, x...)
    end 
end


f = FooBar() 

println(f[0, 5]) # Returns 100 Correctly
println(f[1, 1]) # Causes stackoverflow, want to return 2