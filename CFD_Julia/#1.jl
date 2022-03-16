# This is the first Julia code on CFD

# Here we seek to solve the 1D Linear Convection where the equation is:

#====================================#
#                                    #
#       ∂u/∂t + c. ∂u/∂x = 0         #
#                                    #
#====================================#

using Plots

# Define the spatial grid

nx = 41 # number of x nodes
dx = 2/(nx-1) # number of each x node size
nt = 25 # number of timesteps
dt = 0.025 # Δt
c = 1 # wavespeed

# Define the initial condition for this PDE --> a hat function

u = ones(nx)

u[Int((0.5/dx)):Int((1/dx+1))] .= 2
x = LinRange(0,2,nx)

display(plot(x,u))

un = ones(nx)

for n in 1:nt
    local un = copy(u)
    for i in 2:nx
        u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])
    end
end

display(plot(x,u))


#= Try to change the scopus of un inside the for loop
consider change it to global and later print the variable in the REPL
to see the difference between the local and the global scopus.
 =#
#= for n in 1:nt
    local un = copy(u)
    for i in 2:nx
        u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])
    end
end =#


    
