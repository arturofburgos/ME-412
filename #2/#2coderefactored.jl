#==================================================================#
# Coding Project 2 - Solving 2D inviscid scalar transport equation #
#==================================================================#

# Equation to be simulated
#====================================#
#                                    #
#  dc/dt + u. dc/dx + v. dc/dy = 0   #
#                                    #
#====================================#

using ProgressMeter
using PyPlot
ion() # Interactive Output
pygui(true) # Plot the PyPlot in Matplotlib window using the REPL in the VSCode, either true or false

# Create the 2D spatial grid

println("Setting up the spatial grid...")
nx = 321 # number of x elements
ny = 321 # number of y elements 
length = 4
dx = length/(nx-1) # element size
dy = length/(ny-1) # element size

x = LinRange(0,length,nx)
y = LinRange(0,length,ny)

c = zeros(nx,ny) # Populate the c matrix with zeros(nx,ny) 


# Parameters for the initial condition of the PDE

x_min = 0.9 / length * nx 
x_max = 1.1 / length * nx
y_min = 0.9 / length * ny 
y_max = 1.1 / length * ny


# Setting up the initial condition

println("Setting up the initial condition...")
for i in 1:nx
    for j in 1:ny
        if x_min < i < x_max && y_min < j < y_max
            c[i,j] = sin(5*pi*(x[i]-0.9))*sin(5*pi*(y[j]-0.9))
        end
    end
end


# Plot the initial configuration 

figinit = figure()
s=surf(x,y,c, cmap="viridis") # using PyPlot
xlabel("x")
ylabel("y")
zlabel("z")
zlim([-0.8,1.3])
title("Initial configuration of the sine wave")
colorbar(s)
display(gcf())


#= figcontour = figure()
contour(x,y,c) =#


# Define the u and v values

u = 0.5
v = 0.5


# Create the time grid

t_final = 4
dt = 0.001
t_domain = dt:dt:t_final


# Initialize save variables to store the values at needed timestep

save = zeros(nx, ny, t_final)


#====================================#
#                                    #
#             Schemes                #
#                                    #
#====================================#


#===========================#
# First Order Upwind Scheme #
#===========================#
# New Progress Method using ProgressMeter
function FOUS(quantity, time_domain, save_var)
    prog = Progress(round(Int,t_final/dt)+1)
    for t in time_domain
        quantityn = copy(quantity)
        for i in 2:nx-1
            for j in 2:ny-1
                quantity_north = quantityn[i,j]
                quantity_east = quantityn[i,j]
                quantity_south = quantityn[i,j-1]
                quantity_west = quantityn[i-1,j]


                # Discretization of the equation

                quantity[i,j] = (quantityn[i,j] - (v*dt*(quantity_north-quantity_south))/dy - (u*dt*(quantity_east-quantity_west))/dx)
            
            end
        end
    
        next!(prog)   

        if t % 1 == 0 # Nice way to define an Integer, therefore we can save at each second of iteration

            save_var[:,:,Int(t)] = copy(quantity) # Note that we do not need to return any value since we are changing the local biding of the save_var array, therefore it will be change in the global scope

        end

    end 
    
end


#============================#
# Second Order Upwind Scheme #
#============================#
# New Progress Method using ProgressMeter
function SOUS(quantity, time_domain, save_var)
    prog = Progress(round(Int,t_final/dt)+1)
    for t in time_domain
        local quantityn = copy(quantity) # YOU CAN DELETE THE LOCAL WORD THIS IS COOL
        for i in 3:nx-1
            for j in 3:ny-1
                quantity_north = (3/2)*quantityn[i,j] - (1/2)*quantityn[i,j-1]
                quantity_east = (3/2)*quantityn[i,j] - (1/2)*quantityn[i-1,j]
                quantity_south = (3/2)*quantityn[i,j-1] - (1/2)*quantityn[i,j-2]
                quantity_west = (3/2)*quantityn[i-1,j] - (1/2)*quantityn[i-2,j]


                # Discretization of the equation

                quantity[i,j] = (quantityn[i,j] - (v*dt*(quantity_north-quantity_south))/dy - (u*dt*(quantity_east-quantity_west))/dx)

            end
        end

         next!(prog)   

        if t % 1 == 0 # Nice way to define an Integer, therefore we can save at each second of iteration

            save_var[:,:,Int(t)] = copy(quantity)

        end

    
    end
end


#==============#
# Quick Scheme #
#==============#
# New Progress Method using ProgressMeter
function Quick(quantity, time_domain, save_var)
    prog = Progress(round(Int,t_final/dt)+1)
    for t in time_domain
        local quantityn = copy(quantity)
        for i in 3:nx-1
            for j in 3:ny-1
                quantity_north = quantityn[i,j] + (3*quantityn[i,j+1]- 2*quantityn[i,j]-quantityn[i,j-1])/8
                quantity_east = quantityn[i,j] + (3*quantityn[i+1,j]- 2*quantityn[i,j]-quantityn[i-1,j])/8
                quantity_south = quantityn[i,j-1] + (3*quantityn[i,j]- 2*quantityn[i,j-1]-quantityn[i,j-2])/8
                quantity_west = quantityn[i-1,j] + (3*quantityn[i,j]- 2*quantityn[i-1,j]-quantityn[i-2,j])/8


                # Discretization of the equation

                quantity[i,j] = (quantityn[i,j] - (v*dt*(quantity_north-quantity_south))/dy - (u*dt*(quantity_east-quantity_west))/dx)

            end
        end

        next!(prog)   

        if t % 1 == 0 # Nice way to define an Integer, therefore we can save at each second of iteration

            save_var[:,:,Int(t)] = copy(quantity)

        end

    end
end
 

#==============#
# UMIST Scheme #
#==============#
# New Progress Method using ProgressMeter
function UMIST(quantity, time_domain, save_var)
    prog = Progress(round(Int,t_final/dt)+1)
    for t in time_domain
        local quantityn = copy(quantity)
        for i in 3:nx-1
            for j in 3:ny-1


                # For North

                if (quantityn[i,j+1]-quantity[i,j]) == 0
                    quantity_north = quantityn[i,j]
                else
                    rn = (quantityn[i,j]-quantityn[i,j-1])/(quantityn[i,j+1]-quantity[i,j])
                    psi_n_1 = 2*rn
                    psi_n_2 = (1+3*rn)/4
                    psi_n_3 = (3+rn)/4
                    psi_n_4 = 2

                    min_n = min(psi_n_1,psi_n_2,psi_n_3,psi_n_4)

                    if min_n < 0
                        min_n = 0
                    end

                    quantity_north = quantityn[i,j] + min_n/2 * (quantity[i,j+1]-quantity[i,j])
                end


                # For East

                if (quantityn[i+1,j]-quantity[i,j]) == 0
                    quantity_east = quantityn[i,j]
                else
                    re = (quantityn[i,j]-quantityn[i-1,j])/(quantityn[i+1,j]-quantity[i,j])
                    psi_e_1 = 2*re
                    psi_e_2 = (1+3*re)/4
                    psi_e_3 = (3+re)/4
                    psi_e_4 = 2

                    min_e = min(psi_e_1,psi_e_2,psi_e_3,psi_e_4)

                    if min_e < 0
                        min_e = 0
                    end
                
                    quantity_east = quantityn[i,j] + min_e/2 * (quantity[i+1,j]-quantity[i,j])
                end


                # For South

                if (quantityn[i,j]-quantity[i,j-1]) == 0
                    quantity_south = quantityn[i,j-1]
                else
                    rs = (quantityn[i,j-1]-quantityn[i,j-2])/(quantityn[i,j]-quantity[i,j-1])
                    psi_s_1 = 2*rs
                    psi_s_2 = (1+3*rs)/4
                    psi_s_3 = (3+rs)/4
                    psi_s_4 = 2

                    min_s = min(psi_s_1,psi_s_2,psi_s_3,psi_s_4)

                    if min_s < 0
                        min_s = 0
                    end

                    quantity_south = quantityn[i,j-1] + min_s/2 * (quantity[i,j]-quantity[i,j-1])
                end


                # For West

                if (quantityn[i,j]-quantity[i-1,j]) == 0

                    quantity_west = quantityn[i-1,j]

                else
                    rw = (quantityn[i-1,j]-quantityn[i-2,j])/(quantityn[i,j]-quantity[i-1,j])
                    psi_w_1 = 2*rw
                    psi_w_2 = (1+3*rw)/4
                    psi_w_3 = (3+rw)/4
                    psi_w_4 = 2


                    min_w = min(psi_w_1,psi_w_2,psi_w_3,psi_w_4)

                    if min_w < 0
                        min_w = 0
                    end

                    quantity_west = quantityn[i-1,j] + min_w/2 * (quantityn[i,j]-quantityn[i-1,j])
                end


                # Discretization of the equation

                quantity[i,j] = (quantityn[i,j] - (v*dt*(quantity_north-quantity_south))/dy - (u*dt*(quantity_east-quantity_west))/dx)

            end
        end

        next!(prog)   

        if t % 1 == 0 # Nice way to define an Integer, therefore we can save at each second of iteration

            save_var[:,:,Int(t)] = copy(quantity)

        end

    end
end


#====================================#
#                                    #
#          Function Call             #
#                                    #
#====================================#

println("Starting to perform the iterations...")

#FOUS(c,t_domain,save)
SOUS(c,t_domain,save)
#Quick(c,t_domain,save)
#UMIST(c,t_domain,save)


#====================================#
#                                    #
#     Plots and Post-processing      #
#                                    #
#====================================#

function plot_figures(save_var)

    for i =1:t_final
        figure()
        surf_plot = surf(x,y,save_var[:,:,i], cmap="viridis")
        xlabel("x")
        ylabel("y")
        zlabel("z")
        title("Configuration at $i seconds of simulation")
        zlim([-0.8,1.3])
        colorbar(surf_plot)
        display(gcf()) # Plot the PyPlot in VSCode window window using the REPL in the VSCode
    end

end

plot_figures(save)