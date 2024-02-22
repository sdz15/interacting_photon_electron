mu_ph = 0;
mu_el = 0;
alpha_ph = .005;
alpha_el = 1;
kappa = 2;
omega = 5;
ang = 0;
time = 2;
N = 200;
mesh = .01;
pos_alice = mu_ph-1;
pos_bob = mu_el+1;
L = pos_alice-mu_el;

times = (0:mesh:time)+1e-6;

initvals = initialvals([alpha_ph,alpha_el],[mu_ph,mu_el],N,2);
yy_photon = Inf(1,N); % arrival time of photon at Alice when photon and electron interact
yy_single_free = Inf(1,N); % arrival time of electron at Alice when photon and electron do NOT interact
yy_single_boundary = Inf(1,N); % arrival time of electron at Alice with the boundary condition in place

ode_num = length(times);

traj_single_free = Inf(N,ode_num);
traj_interacting_photon = Inf(N,ode_num);
traj_single_boundary = Inf(N,ode_num);

for x=1:N
    q0 = initvals(:,x);
    [tt_single_free,qq_single_free] = ode45(@(t,q) velocity_single_electron(abs(phiMinusApprox(t,q,alpha_el,kappa,0,omega)).^2,abs(phiPlusApprox(t,q,alpha_el,kappa,0,omega)).^2),times,q0(2)-mu_el);
    [tt_single_boundary,qq_single_boundary] = ode45(@(t,q) velocity_single_electron(abs(psiMinusApprox(t,q,alpha_el,kappa,0,omega,abs(L))).^2,abs(psiPlusApprox(t,q,alpha_el,kappa,0,omega,abs(L))).^2),times,q0(2)-mu_el);

    % SAVING VALUES OF TRAJECTORIES TO BE PLOTTED LATER
    traj_single_free(x,:) = qq_single_free;
    traj_single_boundary(x,:) = qq_single_boundary;
    
    f_single_free = find(traj_single_free(x,:)<=L,1);
    f_single_boundary = find(traj_single_boundary(x,:)<=L,1);

    if (~(f_single_free == Inf))
        yy_single_free(x) = tt_single_free(f_single_free);
    end
    if (~(f_single_boundary == Inf))
        yy_single_boundary(x) = tt_single_boundary(f_single_boundary);
    end
end

% COMPUTING THEORETICAL STATISTICS
mesh = .01;
step = 1/mesh;
times_prob = (0:mesh:time)+1e-6;
space = (pos_alice:mesh:pos_bob);

mu_fun_single_free = zeros(1,time/mesh+1);
mu_fun_single_boundary = zeros(1,time/mesh+1);

for t = 1:time/mesh+1
    mu_fun_single_free(t) = -2*j1(abs(phiMinusApprox(times_prob(t),L,alpha_el,kappa,0,omega)).^2,abs(phiPlusApprox(times_prob(t),L,alpha_el,kappa,0,omega)).^2);
    % mu_fun_single_boundary(t) = abs(psiPlus(times_prob(t),L,theta_el,alpha_el,k_el,0,omega,abs(L),step)).^2;
end

txt = {strcat('sigma=',string(alpha_el)),strcat('kappa=',string(kappa)),strcat('omega=',string(omega))};

figure(1)
hold on;
% xlim([pos_alice pos_bob])
xlim([-time time])
ylim([0 time])

for x=1:N
    plot(traj_single_free(x,:),tt_single_free,'Color',[0.8500 0.3250 0.0980]);
    drawnow
    plot(traj_single_boundary(x,:),tt_single_boundary,'Color',[0 0.4470 0.7410]);
    drawnow
end

xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories');
legend({'noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off

% figure(2)
% hold on;
% scatter(initvals(2,:),yy_single_free);
% xlabel('Initial position','FontSize',20);
% ylabel('Arrival time','FontSize',20);
% title('Arrival time of photon vs electron at Alice');
% xlim([pos_alice,pos_bob])
% ylim([0 time])
% hold off
% 
% figure(3)
% hold on;
% cdfplot(yy_single_free);
% plot(times_prob,cumtrapz(times_prob,mu_fun_single_free));
% xlabel('Time','Fontsize',20)
% title('Cumulative distributions of arrival times at Alice');
% xlim([0 time])
% ylim([0 2])
% hold off