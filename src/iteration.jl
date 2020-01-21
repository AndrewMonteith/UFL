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

export pre_order_traversal, pre_order_traversal_, post_order_traversal

struct PreOrderTraversal 
    lifo::Array{AbstractExpr}
end

function Base.iterate(pret::PreOrderTraversal)
    e = pop!(pret.lifo)
    append!(pret.lifo, ufl_operands(e))
    (e, nothing)
end

function Base.iterate(pret::PreOrderTraversal, _::Nothing)
    isempty(pret.lifo) && return nothing 

    Base.iterate(pret)
end

function pre_order_traversal(e::AbstractExpr)
    PreOrderTraversal([e])
end

function pre_order_traversal_(func, op::AbstractExpr)
    to_visit::Array{AbstractExpr} = [op] 

    while !isempty(to_visit)
        e = pop!(to_visit)
        append!(to_visit, ufl_operands(e))
        func(e)
    end
end

struct PostOrderTraversal 
    lifo::Array{Tuple{AbstractExpr, Array{AbstractExpr}}}
end 

function recurse_into_tree!(post::PostOrderTraversal)
    (expr, deps) = post.lifo[end]

    if isempty(deps)
        pop!(post.lifo)
        return expr 
    end 

    while true
        (_, deps) = post.lifo[end]

        dep = pop!(deps)
        dep_ops = ufl_operands(dep)

        if dep_ops === () 
            return dep 
        else 
            push!(post.lifo, (dep, collect(dep_ops)))
        end
    end
end

function Base.iterate(post::PostOrderTraversal)
    leaf_node = recurse_into_tree!(post)
    (leaf_node, nothing)
end

function Base.iterate(post::PostOrderTraversal, _::Nothing)
    isempty(post.lifo) && return nothing 

    Base.iterate(post)
end

function post_order_traversal(e)
    #=
        post order traversal as defined in traversal.py 
        Does not use reversed operands as we don't sort them 
    =#
    PostOrderTraversal([(e, (collect ∘ ufl_operands)(e))])
end