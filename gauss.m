function val = gauss(s,sigma,mu,k)
    val = ((1/(2*pi*sigma))^(1/4)).*(exp(-(((s-mu).^2)/(4*sigma)))).*(cos(k.*s)+(1i*sin(k*s)));
end