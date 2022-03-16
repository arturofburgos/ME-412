# This is the first Julia code on CFD

# Here we seek to solve the Burguers Equation where the equation is:

#===========================================#
#                                           #
#    ∂u/∂t + u. ∂u/∂x = v. ∂^2(u)/∂(x)^2    #
#                                           #
#===========================================#

# Here we will use Symbolics 

using Symbolics
using Plots
using StaticArrays
using LinearAlgebra

@variables x, nu, t # Defining here the symbolic variables
phi = exp(-(x-4*t)^2/(4*nu*(t+1))) + exp(-(x-2*pi-4*t)^2/(4*nu*(t+1)))

D = Differential(x) # Defining the differential operator on x
D(phi) # Applying the operator on phi expression
phiprime = expand_derivatives(D(phi)) # Getting the proper result on the operation, try print this line
u = -2*nu*(phiprime/phi) + 4 # Pluging the previous result and plugging it into the u boundary conditon

ufunc = build_function(u,[t,x,nu]) # Making u a fucntion with t,x and nu as arguments

myf = eval(ufunc)

out = myf(SA[1,4,3])

nx = 101
nt = 100
dx = 2*pi/(nx-1)

nu = 0.07
x = LinRange(0,2*pi,nx)
un = Array{Float64,nx}
t = 0

# Issue from here: make ufunc able to be used into an array not just a single value
# un = [ufunc(t,x0,nu) for x0 in x]

ufunction = build_function(u,[t,x,nu], parallel=Symbolics.MultithreadedForm())
myfunc = eval(ufunction)
outt = myfunc(SA[t, x, nu])

# Try this above with Nick

