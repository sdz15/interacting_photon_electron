function val = gauss(s,alpha,mu,k)
    val = ((1/(2*pi*alpha))^(1/4)).*(exp(-(((s-mu).^2)/(4*alpha)))).*(cos(k.*s)+(1i*sin(k*s)));
end