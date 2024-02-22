function val = phiPlus(t,s,theta,alpha,k,mu,omega,step)
    % X = s-t:((s+t)-(s-t))/step:s+t;

    x1 = psiPlusR(s+t,theta,alpha,k,mu);
    % x2 = -(omega/2)*trapz(X,besselfun(radical(t,s,X),omega,1).*(t-s+X).*psiPlusR(X,theta,alpha,k,mu));
    x2 = -(omega/2)*integral(@(sigma) besselfun(radical(t,s,sigma),omega,1).*(t-s+sigma).*psiPlusR(sigma,theta,alpha,k,mu),s-t,s+t);
    % x3 = -(1i*omega/2)*trapz(X,besselfun(radical(t,s,X),omega,0).*psiMinusR(X,theta,alpha,k,mu));
    x3 = -(1i*omega/2)*integral(@(sigma) besselfun(radical(t,s,sigma),omega,0).*psiMinusR(sigma,theta,alpha,k,mu),s-t,s+t);
    
    val = x1+x2+x3;
end