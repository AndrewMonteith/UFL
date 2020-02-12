macro dump_me(expr)
    sig_exprs = []
    
    functions = findall(expr -> expr isa Expr && expr.head === :function, expr.args[3].args)
    for func_i in functions 
        inner_ctor = expr.args[3].args[func_i]

        # Each Inner Constructor ends with a a call to new(parameter1, ..., parametern)
        new_call_params = inner_ctor.args[2].args[end].args
        
        sig_param_is = findall(param -> param isa Expr && param.args[1] === Symbol("@sig"), new_call_params)
        for sig_param_i ∈ sig_param_is
            # sig_param_i will be the index of parameter wrapped in a @sig 
            # we unwrap it and mark that expression needs to be included in the hash
            new_call_params[sig_param_i] = new_call_params[sig_param_i].args[3]
            push!(sig_exprs, new_call_params[sig_param_i])
        end
    end

    hash_expr = :(hash_me($(sig_exprs...)))

    println(hash_expr)

    expr
end

println(@macroexpand @dump_me struct X
    x::Int 
    y::Bool

    function X()
        println("helllo there?")
        new(@sig(1), @sig(foobar()))
    end

    function X()
        new(1, 2)
    end
end)