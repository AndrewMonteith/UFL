export SpatialCoordinate, find_geometric_dimension

@ufl_type struct SpatialCoordinate <: Terminal 
    ufl_fields=(shape,)

    geometric_dimesnion::Dimension 
    topological_dimension::Dimension 

    function SpatialCoordinate(gd::Dimension, td::Dimension)
        new((gd,), @sig(gd), @sig(td))
    end
end 
Base.show(io::IO, sc::SpatialCoordinate) = print(io, "x")

geometric_dimension(s::SpatialCoordinate) = s.geometric_dimesnion


function find_geometric_dimension(x::AbstractExpr)
    s = Set{Dimension}()

    found_dimension = -1

    @UFL.pre_order_traversal for node ∈ x 
        d = geometric_dimension(node)

        if d === -1 
            continue 
        elseif found_dimension === -1 
            found_dimension = d
        elseif d != found_dimension
            error("Cannot determine geometric dimension from expression")
        end
    end

    found_dimension
end