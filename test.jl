struct Foobar 
    x::String 
end 

Base.:(==)(f::Foobar, f2::Foobar) = f.x[1] === f2.x[1]

f1 = Foobar("ha")
f2 = Foobar("ht")

println(f1 == f2)