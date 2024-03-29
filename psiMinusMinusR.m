function val = psiMinusMinusR(s_ph,s_el,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el)
    val = cos(theta_ph).*gauss(s_ph,alpha_ph,mu_ph,k_ph) .* cos(theta_el).*gauss(s_el,alpha_el,mu_el,k_el);
end