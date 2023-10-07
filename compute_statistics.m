theta_ph = 0;
theta_el = pi/2;
mu_ph = 0;
mu_el = 1;
sigma_ph = .1;
sigma_el = .1;
k_ph = 0;
k_el = 0;
omega = 2;
ang = 0;
time = 4;
N = 50;
bound = 1;
mesh = .01;
step = 1/mesh*10;
pos_alice = mu_ph-1;
pos_bob = mu_el+1;

times = (1e-1:mesh:time)+1e-6;

initvals = initialvals([sigma_ph,sigma_el],[mu_ph,mu_el],N,2);
yy_photon_alice = Inf(1,N); % arrival time of photon at Alice when photon and electron interact
yy_electron_alice = Inf(1,N); % arrival time of electron at Alice when photon and electron do NOT interact

ode_num = time*100-9; %this number may change. Need to see how this number is calculated

traj_noninteracting_alice = Inf(N,ode_num);
traj_noninteracting_bob = Inf(N,ode_num);
traj_interacting_alice = Inf(N,ode_num);
traj_interacting_bob = Inf(N,ode_num);

cmap = colormap;

for x=1:N
    x
    % opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    q0 = initvals(:,x);
    [tt_noninteracting,qq_noninteracting] = ode45(@(t,q) velocity(abs(psiMinusMinus(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2,abs(psiMinusPlus(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2,abs(psiPlusMinusFar(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2,abs(psiPlusPlusFar(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2),times,q0);
    [tt_interacting,qq_interacting] = ode45(@(t,q) velocity(abs(psiMinusMinus(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2,abs(psiMinusPlus(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2,abs(psiPlusMinusNear(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2,abs(psiPlusPlusNear(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2),times,q0);


    % SAVING VALUES OF TRAJECTORIES TO BE PLOTTED LATER
    traj_noninteracting_alice(x,:) = qq_noninteracting(:,1);
    traj_noninteracting_bob(x,:) = qq_noninteracting(:,2);
    traj_interacting_alice(x,:) = qq_interacting(:,1);
    traj_interacting_bob(x,:) = qq_interacting(:,2);

    f_photon = find(traj_interacting_alice(x,:)<=pos_alice,1);
    f_electron_alice = find(traj_noninteracting_bob(x,:)<=pos_alice,1);

    if (~(f_photon == Inf))
        yy_photon_alice(x) = tt_interacting(f_photon);
    end
    if (~(f_electron_alice == Inf))
        yy_electron_alice(x) = tt_noninteracting(f_electron_alice);
    end
end

% % COMPUTING THEORETICAL STATISTICS
% space = (pos_alice:mesh:pos_bob);
% times = (0:mesh:time);
% 
% fun_photon_alice = zeros(time/mesh+1,(mu_el-mu_ph)/mesh+1);
% fun_electron_alice = zeros(time/mesh+1,(mu_el-mu_ph)/mesh+1);
% mu_fun_photon_alice = zeros(1,time/mesh+1);
% mu_fun_electron_alice = zeros(1,time/mesh+1);
% 
% for t = 1:time/mesh+1
%     for s = 1:(pos_bob-pos_alice)/mesh+1
%         % fun(t,s) = j10(abs(psiMinusMinus(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiMinusPlus(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiPlusMinusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2,abs(psiPlusPlusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2)+j01(abs(psiMinusMinus(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiMinusPlus(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiPlusMinusNear(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2,abs(psiPlusPlusNear(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2);
% 
%         % PLOTTING ALICE USING FIRST TWO TERMS OF (31) FROM ARRIVAL TIMES PAPER
%         fun_photon_alice(t,s) = abs(psiPlusMinusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2+abs(psiPlusPlusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;
% 
%         % PLOTTING BOB USING SECOND TWO TERMS OF (31) FROM ARRIVAL TIMES PAPER
%         fun_electron_alice(t,s) = abs(psiMinusPlus(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2+abs(psiPlusPlusFar(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
%     end
%     mu_fun_photon_alice(t) = trapz(fun_photon_alice(t,:))*mesh;
%     mu_fun_electron_alice(t) = trapz(fun_electron_alice(t,:))*mesh;
% end

% TO PLOT THESE STATISTICS, RUN THE PLOT_STATISTICS.M FILE