export FunctionSpace, VectorFunctionSpace

#=
Add required data as needed.
=#
@attach_hash_operators struct FunctionSpace 
    ufl_shape::DimensionTuple
    mesh::@opt_t(Mesh)
    element::AbstractFiniteElement
    
    # degree::Dimension Not sure how to set this... 
   
    # mesh is the domain 
    function FunctionSpace(element::AbstractFiniteElement; @opt(mesh::Mesh))
        shape = if isa(element, MixedElement) 
            error("cannot create FunctionSpace from MixedElement")
        elseif isa(element, VectorElement) 
            tuple(fem_value_shape(element)[:1])
        else
            ()
        end 
        
        new(shape, mesh, element)
    end
end
ufl_element(fs::FunctionSpace) = fs.element
hash_data(fs::FunctionSpace) = ("FunctionSpace", hash_data(fs.mesh), hash_data(fs.element))
Base.show(io::IO, fs::FunctionSpace) = print(io, "FunctionSpace($mesh, $element)")

function VectorFunctionSpace(mesh::Mesh, family::String, @opt(degree::Dimension), @opt(dim::Dimension))
    if dim === nothing 
        dim = geometric_dimension(mesh)
    end

    sub_element = FiniteElement(family; cell=mesh.cell, degree=degree)
    element = VectorElement(sub_element; dim=dim)

    FunctionSpace(element; mesh=mesh)
end