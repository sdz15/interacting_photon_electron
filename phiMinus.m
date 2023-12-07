function val = phiMinus(t,s,theta,sigma,k,mu,omega,step)
    X = s-t:((s+t)-(s-t))/step:s+t;

    x1 = psiMinusR(s-t,theta,sigma,k,mu);
    x2 = -(omega/2)*trapz(X,besselfun(radical(t,s,X),omega,1).*(t+s-X).*psiMinusR(X,theta,sigma,k,mu));
    x3 = -(1i*omega/2)*trapz(X,besselfun(radical(t,s,X),omega,0).*psiPlusR(X,theta,sigma,k,mu));
    
    val = x1+x2+x3;
end