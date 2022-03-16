using SymPy

t,x,nu = symbols("t,x,nu", real=true)
phi = exp(-(x-4*t)^2/(4*nu*(t+1))) + exp(-(x-4*t-2*pi)^2/(4*nu*(t+1)))
phiprime = diff(phi,x)
u  = -2*nu*(phiprime/phi) + 4 
ufunc = lambdify(u,(t,x,nu))
println(ufunc(1,4,3))