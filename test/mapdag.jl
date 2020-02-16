using UFL, Test 

i, j, k = Identity(3), Identity(3), Identity(3)

s = ((i + j) + (i + j)) + k
println(s)

mutable struct Traverser <: Function
    num_of_calls::Int
end 

safe_sum(x::Tuple{}) = 0
safe_sum(x) = sum(x)

function (t::Traverser)(x::UFL.AbstractExpr, operands::UFL.VarTuple{Int})
    t.num_of_calls += 1
    1 + safe_sum(operands)
end

t = Traverser(0) 

@test map_expr_dag(s, t, Int) === 9 
@test t.num_of_calls === 4