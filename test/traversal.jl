using UFL, Test 

i = Identity(3)
j = Identity(3)
k = Identity(3)

s1 = i + j 
F = s1 + k

expected_pre_order_walk = [F, k, s1, i, j]

function pre_t_for(tree)
    walk = []
    UFL.@pre_order_traversal for x in tree
        push!(walk, x)
    end
    walk
end
@test pre_t_for(F) == expected_pre_order_walk

function pre_t_inline(tree)
    walk = []
    UFL.@pre_order_traversal(tree, begin 
        push!(walk, e)
    end)
    walk 
end
@test pre_t_inline(F) == expected_pre_order_walk



expected_post_order_walk = [k, j, i, s1, F]

function post_t_for(tree)
    walk = [] 
    UFL.@post_order_traversal(tree, begin 
        push!(walk, e)
    end)
    walk 
end
@test post_t_for(F) == expected_post_order_walk
    
s2 = s1 + s1 
function pre_unqiue_t(tree)
    walk = []
    UFL.unique_pre_traversal(tree, node -> push!(walk, node))
    walk
end

@test pre_unqiue_t(s2) == [s2, s1, i]