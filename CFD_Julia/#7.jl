# This is the first Julia code on CFD

# Here we seek to solve the 1D Linear Convection where the equation is:

#====================================#
#                                    #
#  ∂u/∂t + c. ∂u/∂x + c. ∂u/∂y = 0   #
#                                    #
#====================================#


using PyPlot
using Plots

nx = 81
ny = 81
nt = 100
c = 1 
dx = 2/(nx - 1)
dy = 2/(ny - 1)
sigma = 0.2
dt = sigma*dx 

x = LinRange(0,2,nx)
y = LinRange(0,2,ny)

u = ones(nx,ny)
un = ones(nx,ny)

u[Int((0.5/dx)):Int((1/dx+1)),Int((0.5/dy)):Int((1/dy+1))] .= 2
fig = figure(figsize=(11,7), dpi=100)
surf(x,y,u, cmap="viridis") # Using PyPlot
#surface(x,y,u) # Using Plots


# Using for loops for spatial discretization
#= for n in 1:nt+1
    global un = copy(u)
    for i in 2:nx-1
        for j in 2:ny-1
            u[i,j] = (un[i,j] - (c*dt/dx*(un[i,j]-un[i,j-1])) - (c*dt/dx*(un[i,j]-un[i-1,j])))

            u[1,:] .= 1 # Dot broadcasting is used to set each value in the matrix slice equall to 1.
            u[end,:] .= 1
            u[:,1] .= 1
            u[:,end] .= 1
        end
    end
end =#

# Using Array iteration
for n in nt+1
    local un = copy(u)
    u[2:end,2:end] = (un[2:end,2:end] - (c*dt/dx*(un[2:end,2:end]-un[2:end,1:end-1])) - (c*dt/dy*(un[2:end,2:end]-un[1:end-1,2:end])))
    u[1,:] .= 1
    u[end,:] .= 1
    u[:,1] .= 1
    u[:,end] .= 1
end

fig = figure(figsize=(11,7), dpi=100)
surf(x,y,u ,cmap="viridis")


