c = [1 2 3 4 5 6 7 8 9 10]

save = ones(Int(size(c,2)/2))

j = 1

for i = 1:10

    if i%2 == 0
        
        save[j] = c[i]
        global j = j+1

    end
    
end

save
