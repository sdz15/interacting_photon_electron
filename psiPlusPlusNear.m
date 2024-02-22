function val = psiPlusPlusNear(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)
    t_til = 1/2*(t_el+s_el-t_ph-s_ph);
    s_til = 1/2*(t_el+s_el+t_ph+s_ph);

    x1 = psiPlusPlusFar(t_ph,s_ph,t_til,s_til,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el,omega,step);

    arr = zeros(1,step);
    i = 1;
    for x = t_til:(t_el-t_til)/step:t_el
        % IF STATEMENT IS TO AVOID NAN ON THE EDGE CASE
        if (s_ph+t_ph-1e-10 <= s_el+t_el-(2*x))
            continue
        end
        arr(i) = psiPlusMinusNear(t_ph,s_ph,x,s_el+t_el-x,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step);
        i = i+1;
    end
    x2 = -1i*omega*trapz((t_el-t_til)/step,arr);

    val = x1+x2;
end