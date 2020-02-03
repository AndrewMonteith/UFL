A sample implemendation of UFL in Julia

### Current Goals

* Make sure indexing works **Done**
* Make add, div, sub, mult work **Sufficiently Done**
* Change API for traversal to make it friendly
* Investigate unique traversal of the tree

### POINTS TO TALK ABOUT

* Is it worth exploring more around traits?

Get this code to work:
```julia 
mesh is 2D

    V = VectorFunctionSpace(mesh, "P", 1)

    v = TestFunction(V) <- Shape (2, )

    u = Function(V) <- pe (2, )

    T = Constant((0, -0.5))

    B = Constant((0, -0.25))

    d = u.geometric_dimension() <- d

    I = Identity(d) (2x2 identity)

    F = I + grad(u) <- (2x2)             # Deformation gradient

    C = F.T*F <- (sum_i F.T[i, j] * F[j, k])                   # Right Cauchy-Green tensor

    # Invariants of deformation tensors

    Ic = tr(C)

    J = det(F)

    # Lamé parameters, "quite squishy"

    mu = Constant(6.3)

    lmbda = Constant(10.0)

    # Stored strain energy density (compressible neo-Hookean model)

    psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

    # Total potential energy
    Pi = psi*dx # - dot(T, u)*ds(4) - dot(B, u)*dx

    F = derivative(Pi, u, v)
```

Think about tranformation pass tools. Don't consider the UFL a tree, consider it a DAG. Think about avoiiding visting the same expression twice. Trees aren't large enough that creating a new tree isn't expensive expensive

### Current Questions for my supervisor
Q: Why is the degree property automatically shadowed in FiniteElement

Q: Should I think about adding Tensors?

Q: Well done for using like the most unknown feature of python (else on for loop)

### Questions & Answers

Q: ufl_domain is commented as saying it should be deprecated? I presume I should do the same.
A: The TODO was wrong

Q: _ufl_hash_data is in UFL and seem to be like "get_unqiue_identifier_for_this_thing". Do I need it
Answer I think: Yes, how else i do "have i seen the same expression"
A: I was right for once

Q: There are two definitions of FiniteElement, one in UFL and one in FIAT. I copied the UFL one. Is that right?
A: Two rights in a row


Q: I've mocked the Mesh class keeping only what we discussed last time but it seems a bit simple.
Q: Trying to rip out all the guf from firedrake into Julia is non-trivial. The point of the project is not the VectorFunctionSpace, Function, ... stuff. It's to focus on the DSL part more? In which case having these trivial copies to merely mock the data is ok?
A: It's fine

Q: How do you determine the degree of a FunctionSpace?
A: Comes from the finite element

Q: Do I need TensorElement/TensorFunctionSpace
A: Add if necessary

Q: elementlist.py How much do I need to copy and deal with?
A: Don't need much from it. properties
Q: VectorFunctionSpace from above example does not exist?

Answer: VectorFunctionSpace lives in Firedrake Project.

Q: What do I need from specific mesh stuff (properties that belong to mesh)

Answer: Need geometric_dimesion & topological dimension
geometric_dimension -> Dimension of space the things lives in (2d if on plane, 3d if in space)
topological_dimension -> Dimension of space needed to represent it (2d for triangle, 3d for cube, ...)
is_affine tag -> We do all our calculations on a reference triangle to simplify stuff. When we map to the phyiscal world (ie our original mesh) the tranformation is affine if it's F(x) = Ax + b

Q: What do the magic strings mean in CellElement and stuff

Answer: Address femtable.org, the string referes the to polynomial(FEM space...) that we use for the triangle.
For symbolic we need shape fo the FEM space
shape of FiniteElement -> () Represents a scalar function over the mes 
shape of VectorElement -> Specified Explicity or (x,) or default gemometric_dimension of mesh 
shape of TensorElement -> Specified Explicity

Q: Parameters in the evaluate function. x, mapping, component & index_values? The evaluate function may not be needed for my side of the project?

Answer: need not worry about evaluate for now

Q: Some methods include a domain parameter like the Constant class. Foudn some examples of domain but they are basically german to me. domain.py?

Answer: Domain is mesh. Constants are defined globally over the entire mesh, not piecewise per triangle

Q: Used in algebra.py, do I have to worry about how to sort Expr objects?

Answer: AN ordering for expressions can make things easier, like common subexpression elimination. It can be considered an optimisation so for now need not worry about it

Q: WHat do I need from function space

Answer: Shape (comes from FEM),
degree of polynomial space,

Q: What do I need from FiniteElement

Answer: Mock square & triangle from UFL. They just provide the properties for the Mesh

Q: Difference between Cell & FiniteElement

Answer: Cell is just a triangle, FE is Cell + Space + Shape + Degree + ...

### Sidenotes
Mesh if just a collection of triangles. Making a finite element specifies what kind of function we want to represent a mesh. FunctionSpace "glues" it togehter, ie describes the collection of the cells
