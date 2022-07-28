# Local and Global biding
# Compute volume

# Tip use SHIFT+ENTER to see the current result for each line 

x = 5

function sphere_1(y)
    
    4/3*π*(y^3)
     
end

sphere_1(x)

@show x # Notice that x did not change its value


#===================================================#
#===================================================#

function sphere_2(y)
    
    4/3*π*(y^3)

    return 12 # Here we show that no matter what the funtion does the return will be crucial
    
end


sphere_2(x)


#===================================================#
#===================================================#



function sphere_3(y)

    return  4/3*π*(y^3) # Here we show that no matter what the funtion does the return will be crucial
    
end

# Or

function sphere_4(y)
    
    4/3*π*(y^3)

    # Here we show that no matter what the funtion does the return will be crucial
    
end



function sphere_5(y)
    
    4/3*π*(y^3)

    nothing  # Here we show that no matter what the funtion does the return will be crucial
    
end


sphere_3(x)

sphere_4(x)

x

x = sphere_3(x) # sphere_4(x)

x

sphere_5(x)

x


#===================================================#
#===================================================#

# Using arrays it is also ok: 

x = [1,2,3,4,5]

function add_in_array(y)

    y = y .+ 1

end


add_in_array(x)

x

x = add_in_array(x)




#===================================================#
#===================================================#

# However when dealing with the index we should be careful: 

x = [1,2,3,4,5]

k = [12,13,14,15,16]

function c1(y)

    y[1] = 112
    
end

c1(x)

@show x


function c2(y)
    y[:] = copy(k[:])

    return y
end


c2(x)

@show x # We have already changed the array here


x = c2(x) # This is not mandatory in order to change the array x 


function c3(y) # Equivalent to c2
    @views y[:] = k[:]
end

c3(x)

x

function c4(y) # Equivalent to c2
    y[:] = k[:]
end

c4(x)