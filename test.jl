struct Foo 
end 

function Base.getindex(f::Foo, x...) 
    println(typeof(x))
end

f = Foo() 

f[1, 1:3]

f[1, 1:3, :]