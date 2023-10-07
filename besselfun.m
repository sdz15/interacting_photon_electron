function val = besselfun(x,omega,index)
    if (norm(x)<=1e-6)
        val = repmat((omega^index)/(factorial(index)*2^index),size(x,1),size(x,2));
        return
    end

    %Removing NaNs from the solution to take care of singularities
    arr = besselj(index,omega*x)./(x.^index);
    arr(isnan(arr))=0;

    val = arr;
end