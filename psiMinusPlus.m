function val = psiMinusPlus(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)
    X = s_el-t_el:((s_el+t_el)-(s_el-t_el))/step:s_el+t_el;

    x1 = psiMinusPlusR(s_ph-t_ph,s_el+t_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el);
    x2 = -(omega/2)*trapz(X,besselfun(radical(t_el,s_el,X),omega,1).*(t_el-s_el+X).*psiMinusPlusR(s_ph-t_ph,X,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el));
    x3 = -(1i*omega/2)*trapz(X,besselfun(radical(t_el,s_el,X),omega,0).*psiMinusMinusR(s_ph-t_ph,X,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el));

    val = x1+x2+x3;
end