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
N = 500;
mesh = .1;
step = 1/mesh;
pos_alice = mu_ph-1;
pos_bob = mu_el+1;

times = (1e-1:mesh:time)+1e-6;

initvals = initialvals([sigma_ph,sigma_el],[mu_ph,mu_el],N,2);
yy_photon_alice = Inf(1,N); % arrival time of photon at Alice when photon and electron interact
yy_electron_alice = Inf(1,N); % arrival time of electron at Alice when photon and electron do NOT interact
yy_single = Inf(1,N); % arrival time of electron at Alice with the boundary condition in place

ode_num = length(times);

traj_noninteracting_alice = Inf(N,ode_num);
traj_noninteracting_bob = Inf(N,ode_num);
traj_interacting_alice = Inf(N,ode_num);
traj_interacting_bob = Inf(N,ode_num);
traj_single = Inf(N,ode_num);

for x=1:N
    x
    q0 = initvals(:,x);
    [tt_noninteracting,qq_noninteracting] = ode45(@(t,q) velocity_helper(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step,false),times,q0);
    [tt_interacting,qq_interacting] = ode45(@(t,q) velocity_helper(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step,true),times,q0);
    [tt_single,qq_single] = ode45(@(t,q) velocity_single_electron(t,q,theta_el,sigma_el,k_el,mu_el,omega,abs(pos_alice),step),times,q0(2));

    % SAVING VALUES OF TRAJECTORIES TO BE PLOTTED LATER
    traj_noninteracting_alice(x,:) = qq_noninteracting(:,1);
    traj_noninteracting_bob(x,:) = qq_noninteracting(:,2);
    traj_interacting_alice(x,:) = qq_interacting(:,1);
    traj_interacting_bob(x,:) = qq_interacting(:,2);
    traj_single(x,:) = qq_single;

    f_photon = find(traj_interacting_alice(x,:)<=pos_alice,1);
    f_electron_alice = find(traj_noninteracting_bob(x,:)<=pos_alice,1);
    f_single = find(traj_single(x,:)<=pos_alice,1);

    if (~(f_photon == Inf))
        yy_photon_alice(x) = tt_interacting(f_photon);
    end
    if (~(f_electron_alice == Inf))
        yy_electron_alice(x) = tt_noninteracting(f_electron_alice);
    end
    if (~(f_single == Inf))
        yy_single(x) = tt_single(f_single);
    end
end

yy_photon_alice(yy_photon_alice==Inf)=NaN;
yy_electron_alice(yy_electron_alice==Inf)=NaN;
yy_single(yy_single==Inf)=NaN;
n_yy_photon_alice=normalize(yy_photon_alice);
n_yy_electron_alice=normalize(yy_electron_alice);
n_yy_single=normalize(yy_single);

COMPUTING THEORETICAL STATISTICS
space = (pos_alice:mesh:pos_bob);
times = (0:mesh:time);

fun_photon_alice = zeros(time/mesh+1,(mu_el-mu_ph)/mesh+1);
fun_electron_alice = zeros(time/mesh+1,(mu_el-mu_ph)/mesh+1);
mu_fun_photon_alice = zeros(1,time/mesh+1);
mu_fun_electron_alice = zeros(1,time/mesh+1);
mu_single = zeros(1,time/mesh+1);

for t = 1:time/mesh+1
    t
    for s = 1:(pos_bob-pos_alice)/mesh+1
        % Wave function for photon j10
        pMM1 = abs(psiMinusMinus(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
        pMP1 = abs(psiMinusPlus(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
        pPM1 = abs(psiPlusMinusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;
        pPP1 = abs(psiPlusPlusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;

        % Wave function for photon j01
        pMM2 = abs(psiMinusMinus(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
        pMP2 = abs(psiMinusPlus(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
        pPM2 = abs(psiPlusMinusNear(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;
        pPP2 = abs(psiPlusPlusNear(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;

        % Wave function for electron j10
        pPM3 = abs(psiPlusMinusFar(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
        pPP3 = abs(psiPlusPlusFar(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;

        % Wave function for electron j10
        pPM4 = abs(psiPlusMinusFar(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
        pPP4 = abs(psiPlusPlusFar(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;

        % Calculating probabilities
        fun_photon_alice(t,s) = -j10(pMM1,pMP1,pPM1,pPP1);
        fun_electron_alice(t,s) = -j01(pMM2,pMP2,pPM4,pPP4);
    end
    mu_fun_photon_alice(t) = trapz(fun_photon_alice(t,:))*mesh;
    mu_fun_electron_alice(t) = trapz(fun_electron_alice(t,:))*mesh;
    mu_single(t) = j0(abs(psiMinus(times(t),-pos_alice,theta_el,sigma_el,k_el,mu_el,omega,abs(pos_alice),step)).^2,abs(psiPlus(times(t),-pos_alice,theta_el,sigma_el,k_el,mu_el,omega,abs(pos_alice),step)).^2)+j0(abs(psiMinus(times(t),pos_alice,theta_el,sigma_el,k_el,mu_el,omega,abs(pos_alice),step)).^2,abs(psiPlus(times(t),pos_alice,theta_el,sigma_el,k_el,mu_el,omega,abs(pos_alice),step)).^2);
end

% TO PLOT THESE STATISTICS, GO TO PLOT_STATISTICS.M