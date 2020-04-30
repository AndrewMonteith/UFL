export ln, sin

abstract type MathFunction <: Operator end 

Base.show(io::IO, m::MathFunction) = print(io, name(m), "(", m.ufl_operands[1], ")")

function count_nodes(root::AbstractExpr)
    i = 0
    @UFL.pre_order_traversal for x ∈ root 
        i += 1 
    end 
    return i
end

@ufl_type struct Ln <: MathFunction 
    ufl_fields = (operands,)
    
    Ln(arg::Real) = ScalarValue(log(arg))

    function Ln(arg::AbstractExpr)
        new(@sig((arg,)))
    end
end 
name(::Ln) = "ln"

ln(x::Real) = Ln(x)
ln(x::AbstractExpr) = Ln(x)

@ufl_type struct Sin <: MathFunction 
    ufl_fields = (operands,)

    function Sin(arg::AbstractExpr)
        !is_true_scalar(arg) && error("sin function was expecting scalar argument")

        new(@sig((arg,)))
    end
end
name(::Sin) = "sin"

sin(x::AbstractExpr) = Sin(x)

@ufl_type struct Cos <: MathFunction 
    ufl_fields = (operands,)

    function Cos(arg::AbstractExpr)
        !is_true_scalar(arg) && error("cos function was expecting scalar argument")

        new(@sig((arg,)))
    end
end
name(::Cos) = "cos"

cos(x::AbstractExpr) = Cos(x)