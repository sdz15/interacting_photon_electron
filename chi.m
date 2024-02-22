function val = chi(t,s,theta,alpha,k,mu,omega,L,step,sign,n)

if t>=(2*n+1)*L-(sign*s)
    X = ((2*n+1)*L-(sign*s)):(t-((2*n+1)*L-(sign*s)))/step:t;

    if sign*((-1)^n)>0
        a1 = trapz(X,psiPlusR(L+t-X,theta,alpha,k,mu).*R(0,sign*s,X,omega,L,n));
        a2 = -(omega/2)*integral2(@(xi,sig) besselfun(radical(t-xi,L,sig),omega,1).*(t-xi-L+sig).*psiPlusR(sig,theta,alpha,k,mu).*R(0,sign*s,xi,omega,L,n),((2*n+1)*L-(sign*s)),t,@(xi) L-t+xi,@(xi) L+t-xi);
        a3 = -(1i*omega/2)*integral2(@(xi,sigma) besselfun(radical(t-xi,L,sigma),omega,0).*psiMinusR(sigma,theta,alpha,k,mu).*R(0,sign*s,xi,omega,L,n),((2*n+1)*L-(sign*s)),t,@(xi) L-t+xi,@(xi) L+t-xi);
        a = -(a1+a2+a3);
    else
        a1 = trapz(X,psiMinusR(-L-t+X,theta,alpha,k,mu).*R(0,sign*s,X,omega,L,n));
        a2 = -(omega/2)*integral2(@(xi,sigma) besselfun(radical(t-xi,-L,sigma),omega,1).*(t-xi-L-sigma).*psiMinusR(sigma,theta,alpha,k,mu).*R(0,sign*s,xi,omega,L,n),((2*n+1)*L-(sign*s)),t,@(xi) -L-t+xi,@(xi) -L+t-xi);
        a3 = -(1i*omega/2)*integral2(@(xi,sigma) besselfun(radical(t-xi,-L,sigma),omega,0).*psiPlusR(sigma,theta,alpha,k,mu).*R(0,sign*s,xi,omega,L,n),((2*n+1)*L-(sign*s)),t,@(xi) -L-t+xi,@(xi) -L+t-xi);
        a = -(a1+a2+a3);
    end
    
    x1 = a;
else 
    x1 = 0;
end

if t>=(2*n+1)*L+(sign*s)
    X = ((2*n+1)*L+(sign*s)):(t-((2*n+1)*L+(sign*s)))/step:t;

    if sign*((-1)^(n+1))>0
        b1 = trapz(X,psiPlusR(L+t-X,theta,alpha,k,mu).*R(1,-sign*s,X,omega,L,n));
        b2 = -(omega/2)*integral2(@(xi,sigma) besselfun(radical(t-xi,L,sigma),omega,1).*(t-xi-L+sigma).*psiPlusR(sigma,theta,alpha,k,mu).*R(1,-sign*s,xi,omega,L,n),((2*n+1)*L+(sign*s)),t,@(xi) L-t+xi,@(xi) L+t-xi);
        b3 = -(1i*omega/2)*integral2(@(xi,sigma) besselfun(radical(t-xi,L,sigma),omega,0).*psiMinusR(sigma,theta,alpha,k,mu).*R(1,-sign*s,xi,omega,L,n),((2*n+1)*L+(sign*s)),t,@(xi) L-t+xi,@(xi) L+t-xi);
        b = -(b1+b2+b3);
    else
        b1 = trapz(X,psiMinusR(-L-t+X,theta,alpha,k,mu).*R(1,-sign*s,X,omega,L,n));
        b2 = -(omega/2)*integral2(@(xi,sigma) besselfun(radical(t-xi,-L,sigma),omega,1).*(t-xi-L-sigma).*psiMinusR(sigma,theta,alpha,k,mu).*R(1,-sign*s,xi,omega,L,n),((2*n+1)*L+(sign*s)),t,@(xi) -L-t+xi,@(xi) -L+t-xi);
        b3 = -(1i*omega/2)*integral2(@(xi,sigma) besselfun(radical(t-xi,-L,sigma),omega,0).*psiPlusR(sigma,theta,alpha,k,mu).*R(1,-sign*s,xi,omega,L,n),((2*n+1)*L+(sign*s)),t,@(xi) -L-t+xi,@(xi) -L+t-xi);
        b = -(b1+b2+b3);
    end
    
    x2 = b*1i;
else 
    x2 = 0;
end

val = (1i^n)*(x1-x2);