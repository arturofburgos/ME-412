# Array operations:

u = [0 1 2 3 4 5]

for i in 2:length(u)
    println(u[i]-u[i-1])
end
# We can perform the same operation using:

u[2:end]-u[1:end-1]




# Let's revisit 2D linear convection:

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

#= @time begin
    for n in 1:nt+1
        local un = copy(u)
        for i in 2:nx-1
            for j in 2:ny-1
                u[i,j] = (un[i,j] - (c*dt/dx*(un[i,j]-un[i,j-1])) - (c*dt/dx*(un[i,j]-un[i-1,j])))

                u[1,:] .= 1 # Dot broadcasting is used to set each value in the matrix slice equall to 1.
                u[end,:] .= 1
                u[:,1] .= 1
                u[:,end] .= 1
            end
        end
    end
end =#

@time begin
    for n in nt+1
        local un = copy(u)
        u[2:end,2:end] = (un[2:end,2:end] - (c*dt/dx*(un[2:end,2:end]-un[2:end,1:end-1])) - (c*dt/dy*(un[2:end,2:end]-un[1:end-1,2:end])))
        u[1,:] .= 1
        u[end,:] .= 1
        u[:,1] .= 1
        u[:,end] .= 1
    end
end


