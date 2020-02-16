export pre_order_traversal, unique_pre_traversal, post_order_traversal

function decompose_for_loop(for_loop)
    iter_var = for_loop.args[1].args[1]
    root_expression = for_loop.args[1].args[2]
    loop_body = for_loop.args[2]

    root_expression, loop_body, iter_var
end

macro pre_order_traversal(for_loop)
    root, loop_body, iter_var =decompose_for_loop(for_loop)

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

macro post_order_traversal(for_loop)
    root, loop_body, iter_var = decompose_for_loop(for_loop)

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

macro unique_pre_traversal(for_loop)
    root, loop_body, iter_var = decompose_for_loop(for_loop)

    esc(quote 
        to_visit::Vector{UFL.AbstractExpr} = [$root]
        visited = Set{UFL.AbstractExpr}()
        
        while !isempty(to_visit)
            $iter_var = pop!(to_visit)
        
            for op ∈ ufl_operands($iter_var)
                op ∈ visited && continue 
                push!(to_visit, op)
                push!(visited, op)
            end 

            $loop_body
        end
    end)
end