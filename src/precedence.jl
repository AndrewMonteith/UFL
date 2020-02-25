export parstr

precedence(o::Operator) = 0
precedence(o::Sum) = 1
precedence(o::IndexSum) = 2
precedence(o::Union{Product, Division}) = 3
precedence(o::Union{Power}) = 4
precedence(o::Indexed) = 5
precedence(o::Terminal) = 6


function parstr(parent, child)
    p_precedence, c_precedence = precedence(parent), precedence(child)

    if precedence(parent) === 0 
        "(" * string(child) * ")"
    elseif p_precedence > c_precedence
        "(" * string(child) * ")"
    else
        string(child) 
    end 
end