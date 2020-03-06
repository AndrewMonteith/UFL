export ExprList 

@ufl_type struct ExprList <: Operator 
    ufl_fields = (operands,)

    function ExprList(operands...)
        new(@sig(operands))
    end 
end 

Base.length(list::ExprList) = length(list.ufl_operands)
Base.show(io::IO, l::ExprList) = print(io, "ExprList( $(join(map(string, l.ufl_operands), ",")) )")
Base.getindex(list::ExprList, i::Int) = list.ufl_operands[i]
Base.iterate(list::ExprList) = iterate(list.ufl_operands)
Base.iterate(list::ExprList, state) = iterate(list.ufl_operands, state)