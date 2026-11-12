#=
Mapping of geothermal 
=#

"""
y = linearmap_01_to_ab(x)

mapping x ∈ (0,1) => y ∈ (a,b)
"""
linearmap_01_to_ab(x,a,b) = (b-a)*x + a

"""
x = linearmap_ab_to_01(y)

mapping y ∈ (a,b) => x ∈ (0,1)
"""
linearmap_ab_to_01(y,a,b) = (y-a)/(b-a) 

"""
x = tanhmap_R_to_01(y)

mapping y ∈ (-∞,+∞) => x ∈ (0,1)
"""
tanhmap_R_to_01(y) = (tanh(y) + 1)/2 

"""
y = tanhmap_01_to_R(x) 

mapping x ∈ (0,1) => y ∈ (-∞,+∞)
"""
tanhmap_01_to_R(x) = atanh(2*x - 1) 


"""
x = sigmoid_R_to_01(y)

mapping y ∈ (-∞,+∞) => x ∈ (0,1)
"""
sigmoid_R_to_01(y) = 1/(1+exp(-y)) 

"""
y = sigmoid_01_to_R(x) 

mapping x ∈ (0,1) => y ∈ (-∞,+∞)
"""
sigmoid_01_to_R(x) = log(x/(1 - x))

"""
x = sigmoid_R_to_ab(y,a,b)

mapping y ∈ (-∞,+∞) => to x ∈ (a,b)
"""
# map_Rtoab(y,a,b) = (tanh(y) + 1)/2 * (b-a) + a 
sigmoid_R_to_ab(y,a,b) = linearmap_01_to_ab(sigmoid_R_to_01(y),a,b)

"""
y= sigmoid_ab_to_R(x,a,b)

mapping from x ∈ (a,b) to y ∈ (-∞,+∞) 
"""
# tanhmap_Rtoab(x,a,b) = atanh(2*(x-a)/(b-a) - 1) 
sigmoid_ab_to_R(x,a,b) = sigmoid_01_to_R(linearmap_ab_to_01(x,a,b))

"""
x = tanhmap_R_to_ab(y,a,b)

mapping y ∈ (-∞,+∞) => to x ∈ (a,b)
"""
# map_Rtoab(y,a,b) = (tanh(y) + 1)/2 * (b-a) + a 
tanhmap_R_to_ab(y,a,b) = linearmap_01_to_ab(tanhmap_R_to_01(y),a,b)

"""
y= tanhmapp_ab_to_R(x,a,b)

mapping from x ∈ (a,b) to y ∈ (-∞,+∞) 
"""
# tanhmap_Rtoab(x,a,b) = atanh(2*(x-a)/(b-a) - 1) 
tanhmap_ab_to_R(x,a,b) = tanhmap_01_to_R(linearmap_ab_to_01(x,a,b))