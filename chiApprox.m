function val = chiApprox(t,s,alpha,kappa,mu,omega,L,sign,n)

if t>=(2*n+1)*L-(sign*s)
    if sign*((-1)^n)>0
        a = integral(@(xi) -phiPlusApprox(t-xi,L,alpha,kappa,mu,omega).*R(0,sign*s,xi,omega,L,n),((2*n+1)*L-(sign*s)),t);
    else
        a = integral(@(xi) -phiMinusApprox(t-xi,-L,alpha,kappa,mu,omega).*R(0,sign*s,xi,omega,L,n),((2*n+1)*L-(sign*s)),t);
    end
    
    x1 = a;
else 
    x1 = 0;
end

if t>=(2*n+1)*L+(sign*s)
    if sign*((-1)^(n+1))>0
        b = integral(@(xi) -phiPlusApprox(t-xi,L,alpha,kappa,mu,omega).*R(1,-sign*s,xi,omega,L,n),((2*n+1)*L+(sign*s)),t);
    else
        b = integral(@(xi) -phiMinusApprox(t-xi,-L,alpha,kappa,mu,omega).*R(1,-sign*s,xi,omega,L,n),((2*n+1)*L+(sign*s)),t);
    end
    
    x2 = b*1i;
else 
    x2 = 0;
end

val = (1i^n)*(x1-x2);