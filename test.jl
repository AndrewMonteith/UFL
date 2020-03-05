struct X 
    a::Int
    b::Int 
end

function reconstruct(x::T; kwargs...) where T
    fields = fieldnames(T)
    new_members = []

    for field ∈ fields
        push!(new_members,  get(kwargs, field, getfield(x, field)))
    end

    T(new_members...)
end 

x = X(1, 3)

println(x)

y = reconstruct(x; a=4, b=5)
println(y)