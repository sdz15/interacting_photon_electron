function val = psiPlusApprox(t,s,alpha,k,mu,omega,L)

if t>=L-s
    f = -phiPlusApprox(t-L+s,L,alpha,k,mu,omega);
else
    f = 0;
end

n = 0;
chisum = 0;

if t > 0
    while (t>(2*n*L))
        chisum = chisum + chiApprox(t,s,alpha,k,mu,omega,L,1,n);
        n = n+1;
    end
end

val = phiPlusApprox(t,s,alpha,k,mu,omega)+f+chisum; 