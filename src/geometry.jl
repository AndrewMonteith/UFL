export SpatialCoordinate, find_geometric_dimension

abstract type AbstractGeometricCellQuantity <: Terminal end 

@ufl_type struct SpatialCoordinate <: AbstractGeometricCellQuantity
    ufl_fields=(shape,)

    geometric_dimesnion::Dimension 
    topological_dimension::Dimension 

    function SpatialCoordinate(gd::Dimension, td::Dimension)
        new((gd,), @sig(gd), @sig(td))
    end
end 
Base.show(io::IO, sc::SpatialCoordinate) = print(io, "x")

geometric_dimension(s::SpatialCoordinate) = s.geometric_dimesnion

@ufl_type struct Jacobian <: AbstractGeometricCellQuantity 
    ufl_fields=(shape,)
    domain::Mesh

    function Jacobian(domain::Mesh)
        sh = (geometric_dimension(domain), topological_dimension(domain))
        new(sh, @sig(domain))
    end 
end
Base.show(io::IO, j::Jacobian) = print(io, "J")

@ufl_type struct JacobianInverse <: AbstractGeometricCellQuantity
    ufl_fields=(shape,)
    domain::Mesh

    function JacobianInverse(domain::Mesh)
        sh = (topological_dimension(domain), geometric_dimension(domain))
        new(sh, @sig(domain))
    end
end
Base.show(io::IO, ji::JacobianInverse) = print(io, "K")

@ufl_type struct JacobianDeterminant <: AbstractGeometricCellQuantity 
    domain::Mesh

    function Jacobian(domain::Mesh)
        new(domain)
    end
end
Base.show(io::IO, jdet::JacobianDeterminant) = print(io, "detJ")


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