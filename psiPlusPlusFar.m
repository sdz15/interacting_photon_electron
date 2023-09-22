function val = psiPlusPlusFar(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)
    x1 = psiPlusPlusR(s_ph+t_ph,s_el+t_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el);
    x2 = -(omega/2)*integral(@(sigma) besselfun(radical(t_el,s_el,sigma),omega,1).*(t_el-s_el+sigma).*psiPlusPlusR(s_ph+t_ph,sigma,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el),s_el-t_el,s_el+t_el);
    x3 = -(1i*omega/2)*integral(@(sigma) besselfun(radical(t_el,s_el,sigma),omega,0).*psiPlusMinusR(s_ph+t_ph,sigma,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el),s_el-t_el,s_el+t_el);
    
    val = x1+x2+x3;
end