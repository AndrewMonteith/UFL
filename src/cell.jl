export Cell, vertex, triangle, quadrilateral 

num_cell_entities = Dict(
    "vertex" => (1,),
    "interval" => (2, 1),
    "triangle" => (3, 3, 1),
    "quadrilateral" => (4, 4, 1),
    "tetrahedron" => (4, 6, 4, 1),
    "hexahedron" => (8, 12, 6, 1)
)

cell_name_to_dimensions = Dict(cellname => length(v)-1 for (cellname, v) in num_cell_entities)

@attach_hash_operators struct Cell
    name::String
    topological_dimension::Dimension 
    geometric_dimension::Dimension 

    function Cell(cellname::String, geometric_dimension::Dimension=nothing)
        topological_dimension = length(num_cell_entities[cellname]) - 1

        if geometric_dimension === nothing 
            geometric_dimension = topological_dimension
        end 

        if topological_dimension > geometric_dimension 
            error("topological dimension cannot be greater than geometric dimension")
        end

        new(cellname, topological_dimension, geometric_dimension)
    end
end
hash_data(c::Cell) = (c.geometric_dimension, c.topological_dimension, c.name)

function Base.show(io::IO, cell::Cell) 
    gdim, tdim = cell.geometric_dimension, cell.topological_dimension
    
    s = cell.name
    if gdim > tdim 
        s *= "$(gdim)D"
    end 

    show(s)
end

function Base.repr(cell::Cell) 
    gdim, tdim = cell.geometric_dimension, cell.topological_dimension

    if gdim === tdim && cell.name in cell_name_to_dimensions
        cell.name 
    else 
        "Cell($(cell.name), $(gdim))"
    end
end


vertex = Cell("vertex", 0)
triangle = Cell("triangle", 2)
quadrilateral = Cell("quadrilateral", 2)