function val = psiMinusR(s,theta,sigma,k,mu)
    val = cos(theta)*gauss(s,sigma,mu,k);
end