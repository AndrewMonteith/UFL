using Test, UFL 

i = Identity(3)
j = Identity(3)
k = Identity(3)

s1 = i + j 
F = s1 + k

expected_pre_order_walk = [F, k, s1, i, j]
pre_order_walk = []
for x in pre_order_traversal(F)
    push!(pre_order_walk, x)
end

@test pre_order_walk == expected_pre_order_walk

expected_post_order_walk = [k, j, i, s1, F]
post_order_walk = []
for x in post_order_traversal(F)
    push!(post_order_walk, x)
end

@test post_order_walk == expected_post_order_walk