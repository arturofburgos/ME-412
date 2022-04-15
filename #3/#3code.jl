#==================================================================#
#  Coding Project 3 - Solving 2D steady heat conduction equation   #
#==================================================================#

# Equation to be simulated
#====================================#
#                                    #
#  d^2T/dx^2 + d^2T/dy^2 = S(x,y)    #
#                                    #
#====================================#
# Where
#================================================#
#                                                #
#  S(x,y) = -(k1^2 + k2^2)*sin(k1*x)*sin(k2*y)   #
#                                                #
#================================================#

using LinearAlgebra
using Statistics
using DelimitedFiles
using Plots
pyplot()




# Create the spatial grid

nx = 20 
ny = 20
length = pi

dx = length/(nx-1)
dy = length/(ny-1)

x = 0:dx:length
y = 0:dy:length


# Wave length

k1 = 4
k2 = k1


# Boundary Conditions

temperature_left = 0
temperature_right = 0  
temperature_top = 0
temperature_bottom = 0 

#temperature_left = sin(k1*x[1]).*sin.(k2.*y)
#temperature_right = sin(k1*x[nx]).*sin.(k2.*y)
#temperature_top = sin.(k1*x).*sin(k2.*y[1])
#temperature_bottom = sin.(k1*x).*sin(k2.*y[ny])


# Residual and error

err = 1
residue = 10^-7
iter = 0


# Initialize Temperature

T_new = zeros(nx,ny)


# A's coefficients

Aw = 1/dx^2
Ae = 1/dx^2
As = 1/dy^2
An = 1/dy^2
Ap = 2/dx^2 + 2/dy^2



#====================================#
#                                    #
#             Schemes                #
#                                    #
#====================================#
# Attention, this code was developed for one scheme at a time, therefore if you want to see the result
# for Jacobi one, you will need to comment all the other schemes. Use ALT+SHIFT+A to do it.

println("Starting to perform the iterations...")

# Save variable to store the iterations and residual convergence
save_var = ones(2,3500)


#===========================#
#         Analytic          #
#===========================#

T_analytic = zeros(nx,ny)


for i = 2:nx-1
    for j = 2:ny-1
        T_analytic[i,j] = sin(k1*x[i])*sin(k2*y[j])
    end
end



#===========================#
#          Jacobi           #
#===========================#

scheme = "Jacobi"
while err > residue
    
    T_old = copy(T_new)

    for i = 2:nx-1
        for j = 2:ny-1
            s = -(k1^2 + k2^2)*sin(k1*x[i])*sin(k2*y[j])
            T_new[i,j] = (Ae*T_old[i+1,j] + Aw*T_old[i-1,j] + An*T_old[i,j+1] + As*T_old[i,j-1] - s)/Ap

        end
    end
    #global err = norm(T_new-T_old,1)
    global err = mean(mean(abs.(T_new-T_old))) # L1 norm
    global iter = iter + 1
    save_var[1,Int(iter+1)] = iter
    save_var[2,Int(iter+1)] = err
    
    println("At the current iteration: ", iter, ", the error is: ", err)
  
end 



#===========================#
#       Gauss-Seidel        #
#===========================#

#= scheme = "Gauss-Seidel"
while err > residue
    
    T_old = copy(T_new)

    for i = 2:nx-1
        for j = 2:ny-1
            s = -(k1^2 + k2^2)*sin(k1*x[i])*sin(k2*y[j])
            T_new[i,j] = (Ae*T_old[i+1,j] + Aw*T_new[i-1,j] + An*T_old[i,j+1] + As*T_new[i,j-1] - s)/Ap

        end
    end
    #global err = norm(T_new-T_old,1)
    global err = mean(mean(abs.(T_new-T_old))) # L1 norm
    global iter = iter + 1
    save_var[1,Int(iter+1)] = iter
    save_var[2,Int(iter+1)] = err


    println("At the current iteration: ", iter, ", the error is: ", err)
  
end =# 



#===========================#
#           SOR             #
#===========================#

#= scheme = "SOR"
w = 1.4

while err > residue
    
    T_old = copy(T_new)

    for i = 2:nx-1
        for j = 2:ny-1
            s = -(k1^2 + k2^2)*sin(k1*x[i])*sin(k2*y[j])
            T_new[i,j] = (Ae*T_old[i+1,j] + Aw*T_new[i-1,j] + An*T_old[i,j+1] + As*T_new[i,j-1] - s)/Ap
            T_gauss = copy(T_new)
            T_new[i,j] = T_old[i,j]+w*(T_gauss[i,j]-T_old[i,j])
        end
    end
    #global err = norm(T_new-T_old,1)
    global err = mean(mean(abs.(T_new-T_old))) # L1 norm
    global iter = iter + 1
    save_var[1,Int(iter+1)] = iter
    save_var[2,Int(iter+1)] = err


    println("At the current iteration: ", iter, ", the error is: ", err)
  
end =# 



#====================================#
#                                    #
#     Plots and Post-processing      #
#                                    #
#====================================#

# Save into a file
writedlm("#3/data/$nx x $ny $scheme $k1.csv", save_var, ", ")

# Contour Plots
fig1 = contour(x, y, T_analytic, fill = true)
display(plot(fig1, title = "Temperature Contour Analytic Solution: $nx x $ny grid, k1 = $k1", titlefontsize = 11))

fig2 = contour(x, y, T_new, fill = true)
display(plot(fig2, title = "Temperature Contour: $nx x $ny grid, k = $k1", titlefontsize = 11))

# Plot single convergence of L1 norm of the residual. Use the Python code to plot all the schemes combined
display(plot(save_var[1,1:Int(maximum(iter))],save_var[2,1:Int(maximum(iter))], yaxis = :log))


# L1 norm of the DISCRETIZATION error
discretization_error_value = mean(mean(abs.(T_new-T_analytic)))
println(discretization_error_value)


discretization_error = (T_new - T_analytic)
fig3 = contour(x,y, discretization_error, fill = true)
display(plot(fig3, title = "Discretization Error: $discretization_error_value", titlefontsize = 11))
