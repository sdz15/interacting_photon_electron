function val = phiPlusApprox(t,s,sigma,kappa,mu,omega)
    v = kappa/sqrt(1+kappa^2);
    % val = -1i/2*exp(-1i*(pi/4-omega*kappa*s))*((f(s-v*t,sigma,mu)*(sqrt(1+kappa^2)-kappa)*exp(-1i*omega*t*sqrt(1+kappa^2)))-(f(s+v*t,sigma,mu)*(sqrt(1+kappa^2)+kappa)*exp(1i*omega*t*sqrt(1+kappa^2))));
    val = (exp(1i*omega*kappa*s)/(2*sqrt(1+kappa^2)))*(exp(-1i*omega*t*sqrt(1+kappa^2))*(sqrt(1+kappa^2)-kappa)*f(s-v*t,sigma,mu)+exp(1i*omega*t*sqrt(1+kappa^2))*(sqrt(1+kappa^2)+kappa)*f(s+v*t,sigma,mu));
end

function val = f(s,sigma,mu)
    val = ((1/(2*pi*sigma))^(1/4)).*(exp(-(((s-mu).^2)/(4*sigma))));
end