# This is the first Julia code on CFD

# Here we seek to solve the 1D Linear Convection where the equation is:

#====================================#
#                                    #
#       ∂u/∂t + c. ∂u/∂x = 0         #
#                                    #
#====================================#

using Plots



function linearconv(nx)
    
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
end

# linearconv(41)
# linearconv(65)
# linearconv(85)

function cfl_linearconv(nx)
    
    dx = 2/(nx-1) # number of each x node size
    nt = 25 # number of timesteps
    c = 1 # wavespeed


    sigma = 0.5
    dt = sigma*dx
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
end

cfl_linearconv(85)
