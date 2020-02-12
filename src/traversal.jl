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

expr_hashes = Dict{AbstractExpr, UInt32}()

function _compute_expr_hash(expr::AbstractExpr)
    lifo::Vector{Tuple{AbstractExpr, Vector{AbstractExpr}}} = [(expr, (collect ∘ ufl_operands)(expr))]
    
    while !isempty(lifo)
        (expr, deps) = lifo[end]
        
        if isempty(deps)
            if expr ∉ keys(expr_hashes)
                expr_hashes[expr] = compute_hash(expr)
            end
            pop!(lifo)
        else
            e = pop!(deps)
            e_ops = ufl_operands(e)
            if isempty(e_ops) 
                expr_hashes[e] = compute_hash(expr)
            else
                push!(lifo, (e, collect(e_ops)))
            end
        end
    end
    
    expr_hashes[expr]
end

# Was going to use WeakKeyDict but apparently 'objects of type Identity' cannot be finalized
# Apparently the types need be mutable? Until memory becomes a problem this can suffice being a dict
function compute_expr_hash(expr::AbstractExpr)
    get!(expr_hashes, expr, _compute_expr_hash(expr))
end

function unique_pre_traversal(root::AbstractExpr, func::T) where T <: Function
    visited = Set{AbstractExpr}()
    lifo::Vector{AbstractExpr} = [root]
    push!(visited, root)
    
    while !isempty(lifo)
        e = pop!(lifo)

        func(e)

        for op ∈ ufl_operands(e)
            if op ∉ visited 
                push!(lifo, op)
                push!(visited, op)
            end 
        end 
    end
end