export ln

abstract type MathFunction <: Operator end 

Base.show(io::IO, m::MathFunction) = print(io, name(m), "(", m.ufl_operands[1], ")")

@ufl_type struct Ln <: MathFunction 
    ufl_fields=(operands,)
    
    Ln(arg::Real) = ScalarValue(log(arg))

    function Ln(arg::AbstractExpr)
        new(@sig((arg,)))
    end
end 
name(::Ln) = "ln"

ln(x::Real) = Ln(x)
ln(x::AbstractExpr) = Ln(x)