using UFL, Test 

i, j, k = Identity(3), Identity(3), Identity(3)

s = ((i + j) + (i + j)) + k

struct CounterMapper <: UFL.AbstractMapper 
    base::UFL.BaseMapper{Int}
    CounterMapper() = new(UFL.BaseMapper{Int}())
end 

safe_sum(::Tuple{}) = 0
safe_sum(x) = sum(x)

function (cm::CounterMapper)(x::UFL.AbstractExpr)
    1 + safe_sum(map(op -> cm[op], ufl_operands(x)))
end 

@test map_expr_dag(CounterMapper(), s) === 9 
