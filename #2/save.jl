
c = [1 2 3 4 5]

save = ones(2)



#= for i in 1:10
    c[i] = c[i] - 1
    println(c[i])
    
    if i % 2 == 0.0
        for j in 1:5
            save[j] = c[i]
        end
    end
end =#


#= for j = 1:2
    for i =1:5
        if i % 2 == 0
            save[j] = c[i]
            break
        end 
        
    end
end =#

j = 1

for i = 1:5

    

    if i%2 == 0
        
        save[j] = c[i]
        global j = j+1
    else
        continue
        
    end



end



save