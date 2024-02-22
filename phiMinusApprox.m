function val = phiMinusApprox(t,s,alpha,kappa,mu,omega)
    v = kappa/sqrt(1+kappa^2);
    % val = -1i/2*exp(1i*(pi/4+omega*kappa*s))*((f(s-v*t,alpha,mu)*exp(-1i*omega*t*sqrt(1+kappa^2)))+(f(s+v*t,alpha,mu)*exp(1i*omega*t*sqrt(1+kappa^2))));
    val = (exp(1i*omega*kappa*s)./(2*sqrt(1+kappa^2))).*(exp(-1i*omega*t*sqrt(1+kappa^2)).*f(s-v*t,alpha,mu)-exp(1i*omega*t*sqrt(1+kappa^2)).*f(s+v*t,alpha,mu));

end

function val = f(s,alpha,mu)
    val = ((1/(2*pi*alpha))^(1/4)).*(exp(-(((s-mu).^2)/(4*alpha))));
end
