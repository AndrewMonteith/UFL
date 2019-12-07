abstract type AbstractConstantValue <: Terminal end 

is_cellwise_constant(::AbstractConstantValue) = True 
ufl_domains(::AbstractConstantValue) = ()