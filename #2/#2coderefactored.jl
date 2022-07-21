#==================================================================#
# Coding Project 2 - Solving 2D inviscid scalar transport equation #
#==================================================================#

# Equation to be simulated
#====================================#
#                                    #
#  dc/dt + u. dc/dx + v. dc/dy = 0   #
#                                    #
#====================================#

using PyPlot
ion() # Interactive Output
pygui(true) # Plot the PyPlot in Matplotlib window using the REPL in the VSCode


using ProgressMeter
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
surf(x,y,c, cmap="viridis") # using PyPlot
xlabel("x")
ylabel("y")
zlabel("z")
title("Initial configuration of the sine wave")
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

save = zeros(nx, ny, 4)



#====================================#
#                                    #
#             Schemes                #
#                                    #
#====================================#
# Attention, this code was developed for one scheme at a time, therefore if you want to see the result
# for First Order Upwind, you will need to comment all the other schemes. Use ALT+SHIFT+A to do it.


println("Starting to perform the iterations...")


#===========================#
# First Order Upwind Scheme #
#===========================#

#= for t in t_domain
    local cn = copy(c)
    for i in 2:nx-1
        for j in 2:ny-1
            c_north = cn[i,j]
            c_east = cn[i,j]
            c_south = cn[i,j-1]
            c_west = cn[i-1,j]


            # Discretization of the equation

            c[i,j] = (cn[i,j] - (v*dt*(c_north-c_south))/dy - (u*dt*(c_east-c_west))/dx)
            
        end
    end
    
    if round(t/dt) % 100 == 0
        println("timestep: ", round(t/dt))
    end

    if t == 1
        save[:,:,1] = copy(c)
    end
    
    if t == 2
        save[:,:,2] = copy(c)
    end
    
    if t == 3
        save[:,:,3] = copy(c)
    end
    
    if t == 4
        save[:,:,4] = copy(c)
    end
end  =#


#============================#
# Second Order Upwind Scheme #
#============================#
# New Progress Method
prog = Progress(round(Int,t_final/dt)+1)
for t in t_domain
    local cn = copy(c) # YOU CAN DELETE THE LOCAL WORD THIS IS COOL
    for i in 3:nx-1
        for j in 3:ny-1
            c_north = (3/2)*cn[i,j] - (1/2)*cn[i,j-1]
            c_east = (3/2)*cn[i,j] - (1/2)*cn[i-1,j]
            c_south = (3/2)*cn[i,j-1] - (1/2)*cn[i,j-2]
            c_west = (3/2)*cn[i-1,j] - (1/2)*cn[i-2,j]


            # Discretization of the equation

            c[i,j] = (cn[i,j] - (v*dt*(c_north-c_south))/dy - (u*dt*(c_east-c_west))/dx)
            
        end
    end


    next!(prog)   

    if t == 1
        save[:,:,1] = copy(c)
    end
    
    if t == 2
        save[:,:,2] = copy(c)
    end
    
    if t == 3
        save[:,:,3] = copy(c)
    end
    
    if t == 4
        save[:,:,4] = copy(c)
    end
end 


#==============#
# Quick Scheme #
#==============#

#= for t in t_domain
    local cn = copy(c)
    for i in 3:nx-1
        for j in 3:ny-1
            c_north = cn[i,j] + (3*cn[i,j+1]- 2*cn[i,j]-cn[i,j-1])/8
            c_east = cn[i,j] + (3*cn[i+1,j]- 2*cn[i,j]-cn[i-1,j])/8
            c_south = cn[i,j-1] + (3*cn[i,j]- 2*cn[i,j-1]-cn[i,j-2])/8
            c_west = cn[i-1,j] + (3*cn[i,j]- 2*cn[i-1,j]-cn[i-2,j])/8


            # Discretization of the equation

            c[i,j] = (cn[i,j] - (v*dt*(c_north-c_south))/dy - (u*dt*(c_east-c_west))/dx)
             
        end
    end
    
    if round(t/dt) % 100 == 0
        println("timestep: ", round(t/dt))
    end

    if t == 1
        save[:,:,1] = copy(c)
    end
    
    if t == 2
        save[:,:,2] = copy(c)
    end
    
    if t == 3
        save[:,:,3] = copy(c)
    end
    
    if t == 4
        save[:,:,4] = copy(c)
    end
end  =#


#==============#
# UMIST Scheme #
#==============#

#= for t in t_domain
    local cn = copy(c)
    for i in 3:nx-1
        for j in 3:ny-1


            # For North

            if (cn[i,j+1]-c[i,j]) == 0
                c_north = cn[i,j]
            else
                rn = (cn[i,j]-cn[i,j-1])/(cn[i,j+1]-c[i,j])
                psi_n_1 = 2*rn
                psi_n_2 = (1+3*rn)/4
                psi_n_3 = (3+rn)/4
                psi_n_4 = 2

                min_n = min(psi_n_1,psi_n_2,psi_n_3,psi_n_4)

                if min_n < 0
                    min_n = 0
                end

                c_north = cn[i,j] + min_n/2 * (c[i,j+1]-c[i,j])
            end

            
            # For East

            if (cn[i+1,j]-c[i,j]) == 0
                c_east = cn[i,j]
            else
                re = (cn[i,j]-cn[i-1,j])/(cn[i+1,j]-c[i,j])
                psi_e_1 = 2*re
                psi_e_2 = (1+3*re)/4
                psi_e_3 = (3+re)/4
                psi_e_4 = 2

                min_e = min(psi_e_1,psi_e_2,psi_e_3,psi_e_4)

                if min_e < 0
                    min_e = 0
                end
            
                c_east = cn[i,j] + min_e/2 * (c[i+1,j]-c[i,j])
            end


            # For South
            
            if (cn[i,j]-c[i,j-1]) == 0
                c_south = cn[i,j-1]
            else
                rs = (cn[i,j-1]-cn[i,j-2])/(cn[i,j]-c[i,j-1])
                psi_s_1 = 2*rs
                psi_s_2 = (1+3*rs)/4
                psi_s_3 = (3+rs)/4
                psi_s_4 = 2

                min_s = min(psi_s_1,psi_s_2,psi_s_3,psi_s_4)

                if min_s < 0
                    min_s = 0
                end

                c_south = cn[i,j-1] + min_s/2 * (c[i,j]-c[i,j-1])
            end


            # For West

            if (cn[i,j]-c[i-1,j]) == 0

                c_west = cn[i-1,j]

            else
                rw = (cn[i-1,j]-cn[i-2,j])/(cn[i,j]-c[i-1,j])
                psi_w_1 = 2*rw
                psi_w_2 = (1+3*rw)/4
                psi_w_3 = (3+rw)/4
                psi_w_4 = 2


                min_w = min(psi_w_1,psi_w_2,psi_w_3,psi_w_4)

                if min_w < 0
                    min_w = 0
                end

                c_west = cn[i-1,j] + min_w/2 * (cn[i,j]-cn[i-1,j])
            end


            # Discretization of the equation

            c[i,j] = (cn[i,j] - (v*dt*(c_north-c_south))/dy - (u*dt*(c_east-c_west))/dx)
            
        end
    end
    
    if round(t/dt) % 100 == 0
        println("timestep: ", round(t/dt))
    end

    if t == 1
        save[:,:,1] = copy(c)
    end
    
    if t == 2
        save[:,:,2] = copy(c)
    end
    
    if t == 3
        save[:,:,3] = copy(c)
    end
    
    if t == 4
        save[:,:,4] = copy(c)
    end
end  =#

#====================================#
#                                    #
#     Plots and Post-processing      #
#                                    #
#====================================#


fig1 = figure()
surf(x,y,save[:,:,1], cmap="viridis") # using PyPlot
xlabel("x")
ylabel("y")
zlabel("z")
title("Configuration at 1 second of simulation")
display(gcf()) # Plot the PyPlot in VSCode window window using the REPL in the VSCode

fig2 = figure()
surf(x,y,save[:,:,2], cmap="viridis") # using PyPlot
xlabel("x")
ylabel("y")
zlabel("z")
title("Configuration at 2 seconds of simulation")
display(gcf()) # Plot the PyPlot in VSCode window window using the REPL in the VSCode

fig3 = figure()
surf(x,y,save[:,:,3], cmap="viridis") # using PyPlot
xlabel("x")
ylabel("y")
zlabel("z")
title("Configuration at 3 seconds of simulation")
display(gcf()) # Plot the PyPlot in VSCode window window using the REPL in the VSCode

fig4 = figure()
display(surf(x,y,save[:,:,4], cmap="viridis")) # using PyPlot
xlabel("x")
ylabel("y")
zlabel("z")
title("Configuration at 4 seconds of simulation")
display(gcf()) # Plot the PyPlot in VSCode window window using the REPL in the VSCode



#function save_info()
