theta = pi/2;
mu = 0;
alpha = 1;
kappa = 1;
omega = 20;
k = kappa*omega;
ang = 0;
time = 10;
mesh = 1/100;
step = 1/mesh;

N = 5;

times = (0:mesh:time)+1e-6;

initvals = initialvals(alpha,mu,N,1);
yy_single_free = Inf(1,N); 

ode_num = length(times);

traj_single_free = Inf(N,ode_num);

for x=1:N
    x
    q0 = initvals(x);
    [tt_single_free,qq_single_free] = ode45(@(t,q) velocity_single_electron(abs(phiMinus(t,q,theta,alpha,k,mu,omega,step)).^2,abs(phiPlus(t,q,theta,alpha,k,mu,omega,step)).^2),times,q0-mu_el);
    [tt_single_free_approx,qq_single_free_approx] = ode45(@(t,q) velocity_single_electron(abs(phiMinusApprox(t,q,alpha,kappa,mu,omega)).^2,abs(phiPlusApprox(t,q,alpha,kappa,mu,omega)).^2),times,q0-mu_el);

    v_diff = zeros(N,ode_num);
    s_diff = zeros(N,ode_num);
    for i=1:ode_num
        pM = phiMinus(times(i),qq_single_free(i),theta,alpha,k,mu,omega,step);
        pP = phiPlus(times(i),qq_single_free(i),theta,alpha,k,mu,omega,step);
        pMApprox = phiMinusApprox(times(i),qq_single_free(i),alpha,kappa,mu,omega);
        pPApprox = phiPlusApprox(times(i),qq_single_free(i),alpha,kappa,mu,omega);

        v_exact = j1(pM,pP)./j0(pM,pP);
        v_spa = j1(pMApprox,pPApprox)./j0(pMApprox,pPApprox);
        v_diff(x,i) = v_exact-v_spa;

        s_diff(x,i) = qq_single_free(i)-qq_single_free_approx(i);
    end
end

figure(1)
hold on
for x=1:N
    plot(times,v_diff(x,:));
end
xlabel('Time','FontSize',20);
ylabel('Velocity error','FontSize',20);
title('Error in velocity with respect to time');
hold off

figure(2)
hold on
for x=1:N
    plot(times,s_diff(x,:));
end
xlabel('Time','FontSize',20);
ylabel('Position error','FontSize',20);
title('Error in position with respect to time');
hold off