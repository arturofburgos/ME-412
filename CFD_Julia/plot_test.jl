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

plot(x,u)

un = ones(nx)


for n in 1:nt
    local un = copy(u)
    for i in 2:nx
        u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])
    end
    plot(x,u)
end

plot(x,u)

@userplot oneconvec
@recipe function f(oc::oneconvec)
    x,y,i,t = oc.args
    title --> "1D Convection"
    xaxis --> ("x", (0,2))
    yaxis --> ("Temperature", (0,2.5))
    seriestype --> :line
    linewidth --> 2
    legend --> :none
    dpi --> 300
    c --> cgrad(:thermal)
    line_z --> y
    annotations --> [(0.1, -0.70, "α = " * string(a)), 
                    (0.137, -0.85, "T = " * @sprintf("%1.3f", t[i]) * " s")]
    annotationfontsize --> 12
    x, y
end