export pre_order_traversal, unique_pre_traversal, post_order_traversal

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

macro post_order_traversal(for_loop)
    iter_var = for_loop.args[1].args[1]
    root_expression = for_loop.args[1].args[2]
    loop_body = for_loop.args[2]

    make_post_traversal_pattern(root_expression, loop_body; iter_var=iter_var)
end

function unique_pre_traversal(root::AbstractExpr, func::T) where T <: Function
    to_visit::Vector{UFL.AbstractExpr} = [root]
    visited = Set{AbstractExpr}()

    while !isempty(to_visit)
        x = pop!(to_visit)

        func(x)

        for op ∈ ufl_operands(x)
            op ∈ visited && continue 
            push!(to_visit, op)
            push!(visited, op)
        end 
    end
end