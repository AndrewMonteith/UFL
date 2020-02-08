export Indexed, IndexSum, ComponentTensor

@ufl_type struct Indexed <: Operator
    ufl_fields = (operands,free_indices, index_dimensions)

    function Indexed(expr::AbstractExpr, multiindex::MultiIndex)
        operands = (expr, as_ufl(multiindex))

        shape = ufl_shape(expr)

        if length(shape) !== length(multiindex)
            error("Invalid number of indices $(length(multiindex)) for tensor expression of rank $(length(shape))")
        elseif any((si < di.d for (si, di) ∈ zip(shape, multiindex) if di isa FixedIndex))
            error("Fixed index out of range")
        end

        efi = ufl_free_indices(expr)
        efid = ufl_index_dimensions(expr)

        # If efi or efid are empty collect destroys the type information of the elements of the array 
        # by erasing them to Union{} which causes the push later to fail
        fi_and_d::Vector{Tuple{AbstractIndex, Dimension}} = if isempty(efi) 
            []
        else
            collect(zip(efi, efid))
        end

        for (pos, ind) ∈ enumerate(multiindex)
            if ind isa Index 
                push!(fi_and_d, (ind, shape[pos]))
            end 
        end 

        unique!(fi_and_d)

        fi, fid = if isempty(fi_and_d)
            (), ()
        else
            collect(zip(fi_and_d...))
        end 

        new(operands, fi, fid)
    end
end

@ufl_type struct IndexSum <: Operator 
    ufl_fields = (operands,free_indices, index_dimensions)
    
    dim::Dimension
    
    function IndexSum(summand::AbstractExpr, index::MultiIndex)
        length(index) !== 1 && error("Expecting a single Index but got $(length(index))")
        
        j, = index 
        fi, fid = ufl_free_indices(summand), ufl_index_dimensions(summand)
        pos = findall(i -> i == j, fi)[1]
        
        new_fi = tuple(fi[1:pos-1]..., fi[pos+1:end]...)
        new_fid = tuple(fid[1:pos-1]..., fid[pos+1:end]...)
        
        new((summand, as_ufl(index)), new_fi, new_fid, pos)
    end
end


function create_slice_indices(indexer, shape, fi)
    all_indices::Array{AbstractIndex} = []
    free_indices::Array{AbstractIndex} = []
    repeated_indices::Array{AbstractIndex} = []
    slice_indices::Array{AbstractIndex} = [] 

    for ind ∈ indexer 
        if ind isa Index 
            push!(all_indices, ind)
            (ind ∈ fi || ind ∈ free_indices) && push!(repeated_indices, ind)
            push!(free_indices, ind)
        elseif ind isa FixedIndex 
            ind.d > shape[length(all_indices) + 1] && error("Index out of bounds.")
            push!(all_indices, ind) 
        elseif ind isa Int 
            ind > shape[length(all_indices) + 1] && error("Index out of bounds.")
            push!(all_indices, FixedIndex(ind)) 
        elseif ind isa Colon 
            ind = Index() 
            push!(slice_indices, ind)
            push!(all_indices, ind)
        else
             error("Not expecting $(ind)")
        end
    end 

    length(all_indices) !== length(shape) && error("Component and shape length don't match")
    
    tuple(all_indices...), tuple(slice_indices...), tuple(repeated_indices...)
end

function Base.getindex(e::AbstractExpr, indexer...)
    #=
        Julia does not handle the elipsis slicing
    =#

    if !(indexer isa Tuple)
        indexer = (indexer,)
    end

    shape = ufl_shape(e)
    
    all_indices, slice_indices, repeated_indices = create_slice_indices(indexer, shape, ufl_free_indices(e))

    # If foo[:] === foo
    length(all_indices) === length(slice_indices) && return e 

    a = Indexed(e, all_indices)

    if !isempty(slice_indices)
        a = as_tensor(a, slice_indices)
    end

    for ii ∈ repeated_indices
        a = IndexSum(a, (ii,))
    end 

    e isa Zero && return Zero(ufl_shape(a), ufl_free_indices(a), ufl_index_dimensions(a))

    return a
end