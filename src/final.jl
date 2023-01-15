"""Quick implementation of LP bounds for sphere packing. In order to 
change the hyperparameters one has to change the hyper dictionary
in main file."""


include("./main.jl")


using .sphere_main,Plots

#  D Euclidean space dimension for sphere packing,
# 2N  maximum number of equation. For bounds to converge N~D
#full_list is a dictionary to keep track of roots
#lpSpacking is the function which solves the equations and give back roots and coefficients
D=24
N=60
full_list= Dict{Int64,Matrix{BigFloat}}()
sphere_main.lpSpacking(24, 60,full_list)

