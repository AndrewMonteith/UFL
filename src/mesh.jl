export Mesh, UnitTriangleMesh, UnitSquareMesh

#=
    Meshes are non-trivial in the FEM world.
    I will add data to this structure as I need it.
=#

# MAYBE TODO:
# - @attach_ufl_id
# - @attach_operators_from_hash_data
# - is_affine
# - Information for actual data of mesh (ie triangle or square?)
struct Mesh 
    geometric_dimension::Dimension 
    topological_dimension::Dimension 

    function Mesh(c::Cell)
        new(geometric_dimension(c), topological_dimension(c))
    end

    function Mesh(geometric_dimension::Dimension, topological_dimension::Dimension)
        new(geometric_dimension, topological_dimension)
    end
end


#=
    in the real project these functions take additional data to properly represent the actual 
    geometrical data. However for now I only mock the dimension properties so need not accept that information
=#
UnitTriangleMesh() = Mesh(triangle)
UnitSquareMesh() = Mesh(quadrilateral)