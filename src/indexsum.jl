export IndexSum

struct IndexSum <: Operator 
  
    function IndexSum(summand::AbstractExpr, index::MultiIndex)
        if length(index) != 1 
            error("expecting MultiIndex with only 1 index")
        end

        # This is where one could do Zero Simplification

        i, = index

        new()
    end

end