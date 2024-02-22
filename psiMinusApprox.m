function val = psiMinusApprox(t,s,alpha,kappa,mu,omega,L)

if t>=L+s
    f = -phiMinusApprox(t-L-s,-L,alpha,kappa,mu,omega);
else
    f = 0;
end

n = 0;
chisum = 0;

if t > 0
    while (t>(2*n*L))
        chisum = chisum + chiApprox(t,s,alpha,kappa,mu,omega,L,-1,n);
        n = n+1;
    end
end

val = phiMinusApprox(t,s,alpha,kappa,mu,omega)+f+chisum;