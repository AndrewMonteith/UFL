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

# struct PreOrderTraversal 
#     lifo::Array{AbstractExpr}
# end

# function Base.iterate(pret::PreOrderTraversal)
#     e = pop!(pret.lifo)
#     append!(pret.lifo, ufl_operands(e))
#     (e, nothing)
# end

# function Base.iterate(pret::PreOrderTraversal, _::Nothing)
#     isempty(pret.lifo) && return nothing 

#     Base.iterate(pret)
# end

# function pre_order_traversal(e::AbstractExpr)
#     PreOrderTraversal([e])
# end

struct PreOrderWalker 
    node::AbstractExpr
end

mutable struct CapacityArray{T} 
    buf::Vector{T} 
    len::Int 
    capacity::Int

    function CapacityArray{T}(capacity::Int) where T
        new{T}(Vector{T}(undef, capacity), 0, capacity)
    end
end

@inline function Base.append!(arr::CapacityArray{T}, vals) where T
    n = length(vals)

    if arr.len + n > arr.capacity
        append!(arr.buf, Vector{T}(undef, arr.capacity))
        arr.capacity = round(arr.capacity * 1.5)
    end

    @inbounds for ii in 1:n 
        arr.buf[arr.len + ii] = vals[ii]
    end

    arr.len += n 
end

@inline function Base.pop!(arr::CapacityArray{T}) where T 
    x = arr.buf[arr.len]
    arr.len -= 1
    x 
end

function new_state(x::AbstractExpr)::Vector{AbstractExpr}
    v = Vector{AbstractExpr}()
    push!(v, x)
    v
end

function Base.iterate(w::PreOrderWalker, stack=new_state(w.node))
    isempty(stack) && return nothing 
    
    e = pop!(stack)
    append!(stack, ufl_operands(e))
    (e, stack)
end

Base.eltype(w::PreOrderWalker) = AbstractExpr

@inline function pre_order_traversal(e::AbstractExpr)::PreOrderWalker
    PreOrderWalker(e)
end

@inline function pre_order_traversal_(func::F, op::AbstractExpr) where F <: Function
    to_visit::Array{AbstractExpr} = [op]

    while !isempty(to_visit)
        e = pop!(to_visit)
        append!(to_visit, ufl_operands(e))
        func(e)
    end
end

@inline function pre_order_traversal_2(func::F, op::AbstractExpr) where F <: Function
    to_visit = CapacityArray{AbstractExpr}(1000)
    append!(to_visit, [op])

    while to_visit.len > 0
        e = pop!(to_visit)
        append!(to_visit, ufl_operands(e))
        func(e)
    end
end

macro pre_order_traversal_m(op, code)
    esc(
        quote
            to_visit = CapacityArray{AbstractExpr}(100)
            append!(to_visit, [$op::AbstractExpr])

            while to_visit.len > 0
                e = pop!(to_visit)
                append!(to_visit, ufl_operands(e))

                $code
            end
        end
    )
end

macro pre_order_traversal_m2(op, code)
    esc(
        quote
            to_visit::Vector{AbstractExpr} = [$op] 

            while !isempty(to_visit)
                e = pop!(to_visit)
                append!(to_visit, ufl_operands(e))

                $code
            end
        end
    )
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

        if isempty(dep_ops)
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