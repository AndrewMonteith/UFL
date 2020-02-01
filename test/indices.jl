using UFL, Test 

i = Index() 
j = Index() 

@test j != i 
@test i < j 
@test i <= j 
@test i <= i