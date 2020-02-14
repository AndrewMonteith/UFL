export compute_hash

compute_hash(::Type{T}, datum) where T = hash((ufl_typecode(T), hash(datum)))
compute_hash(::Type{T}, data...) where T = hash((ufl_typecode(T), sort(map(hash, data)...)))
# expr_hashes = Dict{AbstractExpr, UInt32}()

# function _compute_expr_hash(expr::AbstractExpr)
#     lifo::Vector{Tuple{AbstractExpr, Vector{AbstractExpr}}} = [(expr, (collect ∘ ufl_operands)(expr))]
    
#     while !isempty(lifo)
#         (expr, deps) = lifo[end]
        
#         if isempty(deps)
#             if expr ∉ keys(expr_hashes)
#                 expr_hashes[expr] = ufl_compute_hash(expr)
#             end
#             pop!(lifo)
#         else
#             e = pop!(deps)
#             e_ops = ufl_operands(e)
#             if isempty(e_ops) 
#                 expr_hashes[e] = ufl_compute_hash(expr)
#             else
#                 push!(lifo, (e, collect(e_ops)))
#             end
#         end
#     end
    
#     expr_hashes[expr]
# end


# function compute_expr_hash(expr::AbstractExpr)
#     get!(expr_hashes, expr, _compute_expr_hash(expr))
# end

