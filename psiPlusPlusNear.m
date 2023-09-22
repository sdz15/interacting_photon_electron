function val = psiPlusPlusNear(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)
    t_til = 1/2*(t_el+s_el-t_ph-s_ph);
    s_til = 1/2*(t_el+s_el+t_ph+s_ph);
    
    x1 = psiPlusPlusFar(t_ph,s_ph,t_til,s_til,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega);

    N = 10;
    arr = zeros(1,N);
    mesh = (t_el-t_til)/N;
    idx = 1;
    for x = t_til:mesh:t_el
        arr(idx) = psiPlusMinusNear(t_ph,s_ph,x,s_el+t_el-x,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang);
        idx = idx+1;
    end
    x2 = -1i*omega*trapz(mesh,arr);
    
    val = x1+x2;
end