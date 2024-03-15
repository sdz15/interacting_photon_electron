theta = pi/2;
mu = 0;
alpha = 1;
kappa = 2;
omega = 20;
k = kappa*omega;
ang = 0;
time = 30;
mesh = 1/100;
step = 1/mesh;

N = 100;

times = (0:mesh:time)+1e-6;

initvals = initialvals(alpha,mu,N,1);
yy_single_free = Inf(1,N); 

ode_num = length(times);

traj_single_free = Inf(N,ode_num);
v_diff = zeros(N,ode_num);
s_diff = zeros(N,ode_num);

for x=1:N
    x
    q0 = initvals(x);
    [tt_single_free,qq_single_free] = ode45(@(t,q) velocity_single_electron(abs(phiMinus(t,q,theta,alpha,k,mu,omega,step)).^2,abs(phiPlus(t,q,theta,alpha,k,mu,omega,step)).^2),times,q0);
    [tt_single_free_approx,qq_single_free_approx] = ode45(@(t,q) velocity_single_electron(abs(phiMinusApprox(t,q,alpha,kappa,mu,omega)).^2,abs(phiPlusApprox(t,q,alpha,kappa,mu,omega)).^2),times,q0-mu);

    for i=1:ode_num
        pM = abs(phiMinus(times(i),qq_single_free(i),theta,alpha,k,mu,omega,step)).^2;
        pP = abs(phiPlus(times(i),qq_single_free(i),theta,alpha,k,mu,omega,step)).^2;
        pMApprox = abs(phiMinusApprox(times(i),qq_single_free(i),alpha,kappa,mu,omega)).^2;
        pPApprox = abs(phiPlusApprox(times(i),qq_single_free(i),alpha,kappa,mu,omega)).^2;

        v_exact = velocity_single_electron(pM,pP);
        v_spa = velocity_single_electron(pMApprox,pPApprox);
        v_diff(x,i) = v_exact-v_spa;

        s_diff(x,i) = qq_single_free(i)-qq_single_free_approx(i);
    end
end

txt = {strcat('alpha=',string(alpha)),strcat('kappa=',string(kappa)),strcat('omega=',string(omega))};

for x=1:N
    c = [rand,rand,rand];
    figure(1)
    plot(times,v_diff(x,:),'Color',c);
    hold on

    figure(2)
    plot(times,s_diff(x,:),'Color',c);
    hold on
end
figure(1)
xlim([0 time])
xlabel('Time','FontSize',20);
ylabel('Velocity error','FontSize',20);
title('Error in velocity with respect to time');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);

figure(2)
xlim([0 time])
xlabel('Time','FontSize',20);
ylabel('Position error','FontSize',20);
title('Error in position with respect to time');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);

figure(3)
for x=1:N
    scatter(initvals(x),s_diff(x,ode_num));
    hold on
end
xlim([mu-2*alpha mu+2*alpha])
xlabel('Initial position')
ylabel('Error')
title('Error as a function of initial position')
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);

hold off