function vals = initialvals(alpha,mu,N,num_particles)
    arr = [];
    for i=1:num_particles
        arr = [arr; repmat(mu(i),1,N)+randn(1,N)*sqrt(alpha(i))];
    end
    
    vals = arr;

end