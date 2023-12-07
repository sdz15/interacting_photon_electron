function val=velocity_helper(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step,interacting)
    psiMM = abs(psiMinusMinus(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
    psiMP = abs(psiMinusPlus(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;

    if (interacting==true && (s_ph+t_ph-1e-10>s_el-t_el))
        psiPM = abs(psiPlusMinusNear(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;
        psiPP = abs(psiPlusPlusNear(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;
    else
        psiPM = abs(psiPlusMinusFar(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
        psiPP = abs(psiPlusPlusFar(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
    end
    
    val = velocity(psiMM,psiMP,psiPM,psiPP);
end