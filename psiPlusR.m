function val = psiPlusR(s,theta,sigma,k,mu)
    val = sin(theta)*gauss(s,sigma,mu,k);
end