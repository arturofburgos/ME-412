clc 
clear 
close all
format long



ui = 7;  # Inlet velocity
Tin = 25;  # Inlet Temperature
v = 1.6*10^(-5);
X = 1; # Lenght of the plate 
Re = (ui*X)/v; # Reynolds Number


Y = 0.025; # horizontal distance of the plate


# Grid generation

imax = 100;
jmax = 50;


dx = X/imax;
dy = Y/jmax;


# Matrix generation 

U = ones(imax, jmax);
V = zeros(imax, jmax);
x = zeros(1,imax); 
y = zeros(1,jmax);

# Boundary Conditions 

U(:,:) = 1;
U(:,1) = 0;
U(1,:) = ui;
U(:,jmax) = ui;
U(1,1) = 0;


V(:,1) = 0;
V(1,:) = 0;
V(1,1) = 0;

iter = 0

# The main code
for o=1:imax
    for j = 2:jmax-1
         for i=1:imax-1
             U(i+1,j) = U(i,j) + (v*(U(i,j+1)-2*U(i,j)+U(i,j-1))*(dx))/((U(i,j)*dy^2)) - ((U(i,j+1)-U(i,j-1))*(V(i,j)*dx)/(2*U(i,j)*dy));
             V(i+1,j) = V(i+1,j-1) - ((0.5*dy*(U(i+1,j)-U(i,j)+U(i+1,j-1)-U(i,j-1)))/dx);
         end
    end
    
    iter = iter + 1
end



contourf(U',25)
colorbar
