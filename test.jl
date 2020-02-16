function foo(::Type{T}, x::U) where {T <: Any, U <: Any}
    nums = Vector{T}()
end 

foo(Int, 2)