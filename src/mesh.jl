export Mesh, UnitTriangleMesh, UnitSquareMesh

#=
    Meshes are non-trivial in the FEM world.
    I will add data to this structure as I need it.
=#

mesh_id = 0

# MAYBE TODO:
# - @attach_ufl_id
# - is_affine
# - Information for actual data of mesh (ie triangle or square?)
# - @attach_ddoperators_from_hash_data
@attach_hash_operators struct Mesh 
    id::Int
    geometric_dimension::Dimension 
    topological_dimension::Dimension 

    function Mesh(c::Cell)
        Mesh(geometric_dimension(c), topological_dimension(c))
    end

    function Mesh(geometric_dimension::Dimension, topological_dimension::Dimension)
        new_id = mesh_id
        global mesh_id = mesh_id + 1 
        new(new_id, geometric_dimension, topological_dimension)
    end
end
Base.repr(m::Mesh) = "Mesh$(m.id)"
hash_data(m::Mesh) = repr(m)
geometric_dimension(m::Mesh) = m.geometric_dimension


#=
    in the real project these functions take additional data to properly represent the actual 
    geometrical data. However for now I only mock the dimension properties so need not accept that information
=#
UnitTriangleMesh() = Mesh(triangle)
UnitSquareMesh() = Mesh(quadrilateral)