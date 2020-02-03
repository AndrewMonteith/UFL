export pre_order_traversal, post_order_traversal

function make_pre_traversal_pattern(root, loop_body; iter_var::Symbol=:e)
    esc(quote
        # Fully Qualifying the type as Main.UFL.AbstractExpr erases the type information
        # at macro expansion sites 
        to_visit::Vector{UFL.AbstractExpr} = [$root]
        
        while !isempty(to_visit)
            $iter_var = pop!(to_visit)
            append!(to_visit, ufl_operands($iter_var))
        
            $loop_body
        end
    end)
end

macro pre_order_traversal(root, loop_body)
    make_pre_traversal_pattern(root, loop_body)
end

macro pre_order_traversal(for_loop)
    iter_var = for_loop.args[1].args[1]
    root_expression = for_loop.args[1].args[2]
    loop_body = for_loop.args[2]

    make_pre_traversal_pattern(root_expression, loop_body; iter_var=iter_var)
end

function make_post_traversal_pattern(root, loop_body; iter_var::Symbol=:e)
    esc(quote
        lifo::Vector{Tuple{UFL.AbstractExpr, Array{UFL.AbstractExpr}}} = [($root, (collect ∘ ufl_operands)($root))]

        while !isempty(lifo)
            (expr, deps) = lifo[end]

            $iter_var = if isempty(deps)
                pop!(lifo)
                expr 
            else
                dep = pop!(deps)
                dep_ops = ufl_operands(dep)

                if !isempty(dep_ops)
                    push!(lifo, (dep, collect(dep_ops)))
                    continue
                end

                dep
            end 

            $loop_body
        end
    end)
end

macro post_order_traversal(root, loop_body)
    make_post_traversal_pattern(root, loop_body)
end

# struct PostOrderTraversal 
#     lifo::Array{Tuple{AbstractExpr, Array{AbstractExpr}}}
# end 

# function recurse_into_tree!(post::PostOrderTraversal)
#     (expr, deps) = post.lifo[end]

#     if isempty(deps)
#         pop!(post.lifo)
#         return expr 
#     end 

#     while true
#         (_, deps) = post.lifo[end]

#         dep = pop!(deps)
#         dep_ops = ufl_operands(dep)

#         if isempty(dep_ops)
#             return dep 
#         else 
#             push!(post.lifo, (dep, collect(dep_ops)))
#         end
#     end
# end

# function Base.iterate(post::PostOrderTraversal)
#     leaf_node = recurse_into_tree!(post)
#     (leaf_node, nothing)
# end

# function Base.iterate(post::PostOrderTraversal, _::Nothing)
#     isempty(post.lifo) && return nothing 

#     Base.iterate(post)
# end

# function post_order_traversal(e)
#     #=
#         post order traversal as defined in traversal.py 
#         Does not use reversed operands as we don't sort them 
#     =#
#     PostOrderTraversal([(e, (collect ∘ ufl_operands)(e))])
# end