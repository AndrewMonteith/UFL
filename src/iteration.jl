#=
    for i in x 
        ...
    end 

    Unrolls to:
    next = iterate(x)
    while next !== nothing 
        (i, state) = next 
        ... 
        next = iterate(x, state)
    end
=#

mutable struct UniquePreTraversal 

end