#-------------------------------------------#
# Coding Project 1 - Heat Diffusion on an Egg
#-------------------------------------------#
close all
clear all
clc

# Properties According to "http://sites.poli.usp.br/pqi/lea/docs/aiche2003.pdf"

disp("The egg constant properties are:")
k = 0.379 #[W/m.K]
cp = 2677.7 #[J/kg.K]
rho = 1131.5 #[kg/m^3]


alpha = k/(cp*rho) # [m^2/s]

disp("\n\n")

# Create the spatial grid using the a,b axis information

disp("Spatial grid properties are:")
a = 0.04 #[m]
b = 0.03 #[m]

nx = 41 # number of x elements
ny = 31 # number of y elements
nz = 31 # number of z elements

dx = 2*a/(nx-1) # element size
dy = 2*b/(ny-1) # element size
dz = dy # element size

x = linspace(-a,a,nx);
y = linspace(-b,b,ny);
z = linspace(-b,b,nz);

# Create the n_th temperature array

T = zeros(nx,ny,nz);

# Create the isolve aiming to solve only the internal ellipsoid nodes (stair-step treatment)
 
isolve = zeros(nx,ny,nz);

# Evaluate the nodes and set the boundary conditions
disp("\n\n")
disp("Evaluating the nodes and setting the boundary conditions...")

for i=1:nx
  for j=1:ny
    for k=1:nz
      e(i,j,k) = (x(i)./a)^2 + (y(j)./b)^2 + (z(k)./b)^2; # If e(i,j,k) is how the sllipsoid is defined
      if e(i,j,k)<1
        isolve(i,j,k)= 1; # Stair-step method
        T(i,j,k) = 300;
      else
        isolve(i,j,k) = 0; # Stair-step method
        T(i,j,k) = 350;
      end
    endfor
  endfor
endfor

disp("\t\t\t\t\t\t\t\tdone")

# Create the time grid

disp("\n\n")
disp("Creating the time grid using the stability criteria...")
disp("\t\t\t\t\t\t\t\tdone")

# The stability criteria was used in order to create dt
dt = (0.45/alpha)*(1/((1/(dx^2)+1/(dy^2)+1/(dz^2))))

disp("\n\n")
disp("The initial temperature at the ellipsoid centroid:")
disp(T((nx-1)/2,(ny-1)/2,(nz-1)/2)),disp("K")


Tnew = zeros(nx,ny,nz); # Initialize the n+1_th temperature array

# Initialize both temperature and time counters aiming to save the values for later plots
t_stop = 700;
t = zeros(t_stop,1); # Store time data
temperature = zeros(t_stop,1); # Store temperature data 
K1 = 0; # Store specific time where temperature reach 330 K
K2 = 0; # Store specific time where temperature reach 337 K
K3 = 0; # Store specific time where temperature reach 349.5 K

disp("\n\n")
disp("Main iteration happening:")
for g=1:t_stop-1
  for i=2:nx-1
    for j=2:ny-1
      for k=2:nz-1
        if isolve(i,j,k) == 1 
          Tnew(i,j,k) = T(i,j,k) + dt*alpha*(((T(i+1,j,k)-2*T(i,j,k)+T(i-1,j,k))/dx^2) + ((T(i,j+1,k)-2*T(i,j,k)+T(i,j-1,k))/dy^2) + ((T(i,j,k+1)-2*T(i,j,k)+T(i,j,k-1))/dz^2));
          T(i,j,k) = Tnew(i,j,k);          
          if round(T((nx-1)/2,(ny-1)/2,(nz-1)/2)) < 350 # here change in case need to plot the countour at different temperatures, see the imagesc plots (colorbar)
            continue
          else
            break
          endif
        else
          continue
        endif
      endfor
    endfor
  endfor
  if round(T((nx-1)/2,(ny-1)/2,(nz-1)/2)) < 350
    temperature(g) = T((nx-1)/2,(ny-1)/2,(nz-1)/2);
    t(g+1) = t(g) + dt;
    if round(T((nx-1)/2,(ny-1)/2,(nz-1)/2)) == 330
      K1 = t(g);
    endif
    if round(T((nx-1)/2,(ny-1)/2,(nz-1)/2)) == 337
      K2 = t(g);
    endif
    if round(T((nx-1)/2,(ny-1)/2,(nz-1)/2)) == 349
      K3 = t(g);
    endif
    if mod(g,100)==0
      disp("\n")
      disp("\t\t\t\t\t\t\t\t..."),disp(g)
    else
      continue
    endif
  else 
    break
  endif
endfor
disp("\n")
disp("\t\t\t\t\t\t\t\tdone")

disp("Temperature at the boundary")
disp(T(41,31,31))

disp("Temperature at the centroid") # Remember to change the value at line 99 to get at the disired time the temperature contour
disp(T((nx-1)/2,(ny-1)/2,(nz-1)/2))


t1 = t(temperature<330 & temperature>0);
t2 = t((temperature>330 & temperature<337) & temperature>0);
t3 = t(temperature>337 & temperature>0);
T1 = temperature(temperature<330 & temperature>0);
T2 = temperature((temperature>330 & temperature<337) & temperature>0);
T3 = temperature(temperature>337 & temperature>0);
plot(t1,T1)
hold on
plot(t2,T2)
hold on
plot(t3,T3)
legend('\fontsize{16}300 K to 330 K','\fontsize{16} 330 K to 337 K','\fontsize{16} 337 K to 349.5 K')
title('\fontsize{22}Centroid temperature within time')
xlabel('\fontsize{22} time')
ylabel('\fontsize{22} temperature')

imagesc(T(:,:,(nz-1)/2)) # Remember to change the value at line 99 to get at the disired time the temperature contour
colorbar
title('\fontsize{22}Temperature contour along center z plane')
