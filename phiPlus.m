function val = phiPlus(t,s,theta,sigma,k,mu,omega,step)
    X = s-t:((s+t)-(s-t))/step:s+t;

    x1 = psiPlusR(s+t,theta,sigma,k,mu);
    x2 = -(omega/2)*trapz(X,besselfun(radical(t,s,X),omega,1).*(t-s+sigma).*psiPlusR(X,theta,sigma,k,mu));
    x3 = -(1i*omega/2)*trapz(X,besselfun(radical(t,s,X),omega,0).*psiMinusR(X,theta,sigma,k,mu));
    
    val = x1+x2+x3;
end