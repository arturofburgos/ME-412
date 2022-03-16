# This is the first Julia code on CFD

# Here we seek to solve the Burguers Equation where the equation is:

#===========================================#
#                                           #
#    ∂u/∂t + u. ∂u/∂x = v. ∂^2(u)/∂(x)^2    #
#                                           #
#===========================================#

# Here we will use SymPy 

using SymPy
using Plots
using StaticArrays
using LinearAlgebra



t,x,nu = symbols("t,x,nu", real=true) # Defining here the symbolic variables
phi = exp(-(x-4*t)^2/(4*nu*(t+1))) + exp(-(x-4*t-2*pi)^2/(4*nu*(t+1)))
phi_prime = diff(phi,x) # Differentiating phi
u  = -2*nu*(phi_prime/phi) + 4 
ufunc = lambdify(u,(t,x,nu)) # Making u a fucntion with t,x and nu as arguments
println(ufunc(1,4,3))

# Another way to do that: 
#= 
@vars x nu t;
phi = (exp(-(x - 4 * t)^2 / (4 * nu * (t + 1))) +
       exp(-(x - 4 * t - 2 * pi)^2 / (4 * nu * (t + 1))));
phiprime = phi.diff(x);
u = -2 *nu*(phiprime / phi)+4;
ufunc = lambdify(u,[t, x, nu])
println(ufunc(1, 4, 3)) 
=#



nx = 101
nt = 100
Δx = 2 * pi / (nx - 1)
nu = 0.07
Δt = Δx * nu


x = LinRange(0,2*pi,nx)
un = zeros(nx)
t = 0;

u = [ufunc(t, x0, nu) for x0 in x];


# Plot initial configuration
plot(x,u)
scatter!(x,u)
plot!(xlims = (0, 7), ylims = (0, 10))

for n in 1:nt+1
    local un = copy(u)
    for i in 2:nx-1      
        u[i] = un[i] - un[i] * Δt / Δx *(un[i] - un[i-1]) 
            + nu * Δt / Δx^2 *
            (un[i+1] - 2 * un[i] + un[i-1])     
    end
    u[1] = un[1] - un[1] * Δt / Δx * (un[1] - un[end-2]) 
        + nu * Δt / Δx^2 *
        (un[2] - 2 * un[1] + un[end-1])
    u[end] = u[1]
end        
        
u_analytical = [ufunc(nt * Δt, xi, nu) for xi in x]


plot(x,u)
scatter!(x,u)
plot!(x,u_analytical)
plot!(xlims = (0, 7), ylims = (0, 10))
