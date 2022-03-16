# This is the first Julia code on CFD

# Here we seek to solve the 1D Diffusion Equation where the equation is:

#====================================#
#                                    #
#    ∂u/∂t - v. ∂^2(u)/∂(x)^2 = 0    #
#                                    #
#====================================#

using Plots

# Define the spatial grid

nx = 41 # number of x nodes
dx = 2/(nx-1) # number of each x node size
nt = 25 # number of timesteps
v = 0.3 # kinematic viscosity
sigma = 0.2 # parameter used in determining the cfl condition
dt = sigma * dx^2/v

# Define the initial condition for this PDE --> a hat function

u = ones(nx)

u[Int((0.5/dx)):Int((1/dx+1))] .= 2
x = LinRange(0,2,nx)

display(plot(x,u))

un = ones(nx)

for n in 1:nt
    local un = copy(u)
    for i in 2:nx-1
        u[i] = un[i] + v * dt / dx^2 * (un[i+1] - 2 * un[i] + un[i-1])
    end
end

display(plot(x,u))