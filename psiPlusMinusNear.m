function val = psiPlusMinusNear(t_ph,s_ph,t_el,s_el,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)
    s_til = s_el-(s_ph+t_ph);

    X_b = 0:((t_el+s_til)/2)/step:(t_el+s_til)/2;
    X_c = 0:((t_el-s_til)/2)/step:(t_el-s_til)/2;

    x1 = psiPlusMinusFar(0,s_ph+t_ph,(t_el+s_til)/2,s_ph+t_ph+(t_el+s_til)/2,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el,omega,step);
    x2 = exp(-1i.*ang).*psiMinusPlus((t_el-s_til)/2,s_ph+t_ph-(t_el-s_til)/2,(t_el-s_til)/2,s_ph+t_ph-(t_el-s_til)/2,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el,omega,step);

    y1 = psiPlusMinusFar(0,s_ph+t_ph,0,s_ph+t_ph,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el,omega,step);
    y2 = exp(-1i.*ang).*psiMinusPlus(0,s_ph+t_ph,0,s_ph+t_ph,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el,omega,step);
    x3 = -1/2*(y1+y2)*besselfun(radical(t_el,s_til,0),omega,0);

    y3 = trapz(X_b,besselfun((t_el-s_til).*(t_el+s_til-(2*X_b)),omega,1).*psiPlusMinusR(s_ph+t_ph,s_ph+t_ph,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el));
    y4 = -(omega/2)*integral2(@(sigma,b) besselfun((t_el-s_til).*(t_el+s_til-(2*b)),omega,1).*besselfun(radical(b,s_ph+t_ph+b,sigma),omega,1).*(2*b+s_ph+t_ph-sigma).*psiPlusMinusR(s_ph+t_ph,sigma,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el),0,(t_el+s_til)/2,s_ph+t_ph,@(b) s_ph+t_ph+(2*b));
    y5 = -(1i*omega/2)*integral2(@(sigma,b) besselfun((t_el-s_til).*(t_el+s_til-(2*b)),omega,1).*besselfun(radical(b,s_ph+t_ph+b,sigma),omega,0).*psiPlusPlusR(s_ph+t_ph,sigma,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el),0,(t_el+s_til)/2,s_ph+t_ph,@(b) s_ph+t_ph+(2*b));
    x4 = -omega*(t_el-s_til)*(y3+y4+y5);

    y3 = trapz(X_c,besselfun((t_el+s_til).*(t_el-s_til-(2*X_c)),omega,1).*exp(-1i.*ang).*psiMinusPlusR(s_ph+t_ph-(2*X_c),s_ph+t_ph-(2*X_c),theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el));
    y4 = -(omega/2)*integral2(@(sigma,c) besselfun((t_el+s_til).*(t_el-s_til-(2*c)),omega,1).*exp(-1i.*ang).*besselfun(radical(c,s_ph+t_ph-c,sigma),omega,1).*(2*c-s_ph-t_ph+sigma).*psiMinusPlusR(s_ph+t_ph-(2*c),sigma,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el),0,(t_el-s_til)/2,s_ph+t_ph,@(c) s_ph+t_ph+(2*c));
    y5 = -(1i*omega/2)*integral2(@(sigma,c) besselfun((t_el+s_til).*(t_el-s_til-(2*c)),omega,1).*exp(-1i.*ang).*besselfun(radical(c,s_ph+t_ph-c,sigma),omega,0).*psiMinusMinusR(s_ph+t_ph-(2*c),sigma,theta_ph,theta_el,alpha_ph,alpha_el,k_ph,k_el,mu_ph,mu_el),0,(t_el-s_til)/2,s_ph+t_ph,@(c) s_ph+t_ph+(2*c));
    x5 = -omega*(t_el+s_til)*(y3+y4+y5);

    val = x1+x2+x3+x4+x5;
end