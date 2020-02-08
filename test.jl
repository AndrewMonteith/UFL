# do_something(::Type{String}) = println("got a string")
# do_something(::Type{Int}) = println("got an int")
do_something(::Type{T}) where T <: Union{String, Int} = println("got a string or int")


do_something(String)