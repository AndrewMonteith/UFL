struct Ones 
    ones::Vector{Int}

    function Ones(n::Int)
        new(fill(1, n))
    end
end 

Base.iterate(one::Ones) = (1, 1)
function Base.iterate(one::Ones, state::Int)
    if state === length(one.ones)
        nothing 
    else
        (1, state += 1) 
    end
end

x = [Ones(10)...]

println(x)