theta_ph = 0;
theta_el = pi/2;
mu_ph = 0;
mu_el = 1;
sigma_ph = .005;
sigma_el = .005;
k_ph = 0;
k_el = 0;
omega = 2;
ang = 0;
time = 4;
N = 200;
mesh = .01;
step = 1/mesh;
pos_alice = mu_ph-1;
pos_bob = mu_el+1;
L = pos_alice-mu_el;

times = (0:mesh:time)+1e-6;

initvals = initialvals([sigma_ph,sigma_el],[mu_ph,mu_el],N,2);
yy_photon = Inf(1,N); % arrival time of photon at Alice when photon and electron interact
yy_single_free = Inf(1,N); % arrival time of electron at Alice when photon and electron do NOT interact
yy_single_boundary = Inf(1,N); % arrival time of electron at Alice with the boundary condition in place

ode_num = length(times);

traj_single_free = Inf(N,ode_num);
traj_interacting_photon = Inf(N,ode_num);
traj_single_boundary = Inf(N,ode_num);

for x=1:N
    x
    q0 = initvals(:,x);
    [tt_single_free,qq_single_free] = ode45(@(t,q) velocity_single_electron(abs(phiMinus(t,q,theta_el,sigma_el,k_el,0,omega,step)).^2,abs(phiPlus(t,q,theta_el,sigma_el,k_el,0,omega,step)).^2),times,q0(2)-mu_el);
    [tt_interacting_photon,qq_interacting_photon] = ode45(@(t,q) velocity_helper(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step,true),times,q0);
    [tt_single_boundary,qq_single_boundary] = ode45(@(t,q) velocity_single_electron(abs(psiMinus(t,q,theta_el,sigma_el,k_el,0,omega,abs(L),step)).^2,abs(psiPlus(t,q,theta_el,sigma_el,k_el,0,omega,abs(L),step)).^2),times,q0(2)-mu_el);

    % SAVING VALUES OF TRAJECTORIES TO BE PLOTTED LATER
    traj_single_free(x,:) = qq_single_free;
    traj_interacting_photon(x,:) = qq_interacting_photon(:,1);
    traj_single_boundary(x,:) = qq_single_boundary;

    f_photon = find(traj_interacting_photon(x,:)<=pos_alice,1);
    f_single_free = find(traj_single_free(x,:)<=L,1);
    f_single_boundary = find(traj_single_boundary(x,:)<=L,1);

    if (~(f_photon == Inf))
        yy_photon(x) = tt_interacting_photon(f_photon);
    end
    if (~(f_single_free == Inf))
        yy_single_free(x) = tt_single_free(f_single_free);
    end
    if (~(f_single_boundary == Inf))
        yy_single_boundary(x) = tt_single_boundary(f_single_boundary);
    end
end

traj_single_free = traj_single_free+mu_el;
traj_single_boundary = traj_single_boundary+mu_el;

yy_photon(yy_photon==Inf)=NaN;
yy_single_free(yy_single_free==Inf)=NaN;
yy_single_boundary(yy_single_boundary==Inf)=NaN;
n_yy_photon=normalize(yy_photon);
n_yy_single_free=normalize(yy_single_free);
n_yy_single_boundary=normalize(yy_single_boundary);

% COMPUTING THEORETICAL STATISTICS
space = (pos_alice:mesh:pos_bob);

fun_photon = zeros(time/mesh+1,(mu_el-mu_ph)/mesh+1);
mu_fun_photon = zeros(1,time/mesh+1);
mu_fun_single_free = zeros(1,time/mesh+1);
mu_fun_single_boundary = zeros(1,time/mesh+1);

for t = 1:time/mesh+1
    t
    for s = 1:(pos_bob-pos_alice)/mesh+1
        % Wave function for photon j10
        pMM1 = abs(psiMinusMinus(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
        pMP1 = abs(psiMinusPlus(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
        pPM1 = abs(psiPlusMinusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;
        pPP1 = abs(psiPlusPlusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;

        % Calculating probabilities
        fun_photon(t,s) = -j10(pMM1,pMP1,pPM1,pPP1);
    end
    mu_fun_photon(t) = trapz(fun_photon(t,:))*mesh;
    mu_fun_single_free(t) = -2*j1(abs(phiMinus(times(t),L,theta_el,sigma_el,k_el,0,omega,step)).^2,abs(phiPlus(times(t),L,theta_el,sigma_el,k_el,0,omega,step)).^2);
    mu_fun_single_boundary(t) = abs(psiPlus(times(t),L,theta_el,sigma_el,k_el,0,omega,abs(L),step)).^2;
end

figure
hold on;
plot(times,cumtrapz(times,mu_fun_single_free));
plot(times,cumtrapz(times,mu_fun_single_boundary));
xlabel('Time','Fontsize',20)
title('Cumulative distributions of arrival times at Alice');
xlim([0 time])
ylim([0 1])
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

% TO PLOT THESE STATISTICS, GO TO PLOT_STATISTICS.M