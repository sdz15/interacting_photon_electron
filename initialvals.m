function vals = initialvals(sigma,mu,N,num_particles)
    arr = [];
    for i=1:num_particles
        arr = [arr; repmat(mu(i),1,N)+randn(1,N)*sqrt(sigma(i))];
    end
    
    vals = arr;

end