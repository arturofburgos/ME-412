#==================================================================#
#  Coding Project 4 - Solving 2D Steady Boundary Layer Flow        #
#==================================================================#

# Equations to be simulated
# Continuity
#====================================#
#                                    #
#          du/dx + dv/dy = 0         #
#                                    #
#====================================#
# X Momentum (non-conservative)
#============================================================#
#                                                            #
#  rho*( u*(du/dx) + v*(du/dy) ) = -dp/dx + mu*(d^2u/dy^2)   #
#                                                            #
#============================================================#


using Plots
pyplot()

# Initial code information

xlength = 1 # Lenght of the plate 
ylength = 0.05 # Plate height
Reyn = 100000 # Reynolds number -> boundary layer thickness will be approximately 1.6 cm
U_inlet = 1
V_inlet = 0
rho = 1 # Density
mu = 1e-5 # Viscosity
dpdx = 0 # or -0.1


# Create the spatial grid

nx = 10000
ny = 50 
dy = ylength/ny
dx = xlength/nx

y = LinRange(0,ylength,ny)


# Define boundary conditions

u_old = zeros(ny+1)
v_old = zeros(ny+1)
u_new = zeros(ny+1)
v_new = zeros(ny+1)
u_new[ny+1] = 1

for i = 2:ny+1
    u_old[i] = U_inlet
    v_old[i] = V_inlet
end


# Initialize save variables

save_var_u = zeros(10001,51)
save_var_v = zeros(10001,51)


#====================================#
#                                    #
#       Discretization Scheme        #
#                                    #
#====================================#

x = 0
for ix = 1:nx
    if x < xlength
        global x = x + dx
        
        # Solve for u_new

        for j = 2:ny
            conv = v_old[j]*(u_old[j+1]-u_old[j-1])/(2*dy)
            diff = mu/(rho*dy*dy) * (u_old[j+1] - 2*u_old[j] + u_old[j-1])
            u_new[j] = u_old[j] + dx/u_old[j] * (- conv + diff - dpdx)
        end

        # Solve for v_new

        for j = 2:ny
            v_new[j] = v_new[j-1] - (dy/(2*dx)) * (u_new[j] + u_new[j-1] - u_old[j] - u_old[j-1])
        end

        for j = 1:ny
            u_old[j] = u_new[j]
            v_old[j] = v_new[j]
        end

    end

    save_var_u[ix,:] = u_new[:]
    save_var_v[ix,:] = v_new[:]
    
end


#====================================#
#                                    #
#     Plots and Post-processing      #
#                                    #
#====================================#


# Contour Plot

fig1 = contour(save_var_u', fill = true)
display(plot(fig1, title = "Contour Plot:  $nx x $ny grid", titlefontsize = 11 , xaxis = "nx, correspond to 1 m", yaxis = "δ, in mm"))


# Profiles v-velocity

display(plot(save_var_v[2000, 1:end-1],y, title = "Profiles of v-velocity for dP/dx = $dpdx, at each designated x position", titlefontsize = 11,label = "0.2 m", legend = :topright, xaxis = "v-velocity [m/s]", yaxis = "y [m]" ))
plot!(save_var_v[4000, 1:end-1], y, label = "0.4 m")
plot!(save_var_v[6000, 1:end-1], y, label = "0.6 m")
plot!(save_var_v[10000, 1:end-1], y, label = "1.0 m")


# Profiles u-vector equaly spaced

display(plot(save_var_u[1000, 1:end-1],y, title = "Profiles of u-velocity for dP/dx = $dpdx", label = "0.1 m", legend = :topleft, xaxis = "u-velocity [m/s]", yaxis = "y [m]"))
plot!(save_var_u[2000, 1:end-1], y, label = "0.2 m")
plot!(save_var_u[3000, 1:end-1], y, label = "0.3 m")
plot!(save_var_u[4000, 1:end-1], y, label = "0.4 m")
plot!(save_var_u[5000, 1:end-1], y, label = "0.5 m")
plot!(save_var_u[6000, 1:end-1], y, label = "0.6 m")
plot!(save_var_u[7000, 1:end-1], y, label = "0.7 m")
plot!(save_var_u[8000, 1:end-1], y, label = "0.8 m")
plot!(save_var_u[9000, 1:end-1], y, label = "0.9 m")
plot!(save_var_u[10000, 1:end-1], y, label = "1.0 m")


# First 0.99 value array, not used later, however this could be useful in the future

fv = zeros(10000)
fv[:] = findfirst.(x-> x>0.99,eachrow(save_var_u[1:end-1,:])) # first values
display(plot(fv)) # Before interpolation


# Interpolate between the first value higher than 0.99 and the imediate ealier

BL = zeros(10000) # Array to be stored the interpolation

for i=1:nx
    for j=2:ny
        if save_var_u[i,j]>=0.99
            if save_var_u[i,j]>0.99

                local a = save_var_u[i,j]-0.99
                local b = 0.99-save_var_u[i,j-1]
                BL[i] = (j-2)*dy + b*dy/(a+b)
            
            elseif save_var_u[i,j]==0.99
            
                BL[i] = (j-1)*dy  

            end

            break
             
        end
    end
end


#====================================#
#                                    #
#         Blasius Solution           #
#                                    #
#====================================#
# For dp/dx = 0 only, the comparison between Blasius Solution and Numerical Solution


X = 0:10000
Blasius = (5.0 * sqrt.(X))/(sqrt((U_inlet)*(mu/rho)))/10000 # /10000 to make it in [mm]


# Plot comparison between Blasius Solution and Numerical Interplated Solution

display(plot(X,Blasius,legend=:bottomright,title = "Comparison Numerical and Blasius solutions for dP/dx = 0", label = "Blasius", xaxis = "nx, correspond to 1 m", yaxis = "δ, in mm" ,titlefontsize = 11 ))
display(plot(BL*1000,legend=:bottomright,title = "Boundary Layer Thickness for dP/dx = $dpdx", label = "Numerical", xaxis = "nx, correspond to 1 m", yaxis = "δ, in mm" ))
plot!(BL*1000, label="Numerical") # BL*1000 to make it in [mm]
