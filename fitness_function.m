function fitness = fitness_function(indivi,f,dx,dy,HPBW,SLL,Nulls)
    n = size(indivi,1);
    fitness = ones(n,1);
    parfor  i=1:n
        fitness(i)=funP(indivi(i,:),f,dx,dy,HPBW,SLL,Nulls);
    end
end