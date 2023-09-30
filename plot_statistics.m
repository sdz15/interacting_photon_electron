theta_ph = 0;
theta_el = pi/2;
mu_ph = 0;
mu_el = 1;
sigma_ph = .1;
sigma_el = .1;
k_ph = 25;
k_el = 0;
omega = 2;
ang = 0;
time = 4;
N = 50;
bound = 1;
mesh = .01;
pos_alice = mu_ph-1;k
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
    % opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    q0 = initvals(:,x);
    [tt_noninteracting,qq_noninteracting] = ode45(@(t,q) velocity(abs(psiMinusMinus(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiMinusPlus(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiPlusMinusFar(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiPlusPlusFar(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2),times,q0);
    [tt_interacting,qq_interacting] = ode45(@(t,q) velocity(abs(psiMinusMinus(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiMinusPlus(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiPlusMinusNear(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2,abs(psiPlusPlusNear(t,q(1),t,q(2),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2),times,q0);

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

space = (pos_alice:mesh:pos_bob);
times = (0:mesh:time);

fun_photon_alice = zeros(time/mesh+1,(mu_el-mu_ph)/mesh+1);
fun_electron_alice = zeros(time/mesh+1,(mu_el-mu_ph)/mesh+1);
mu_fun_photon_alice = zeros(1,time/mesh+1);
mu_fun_electron_alice = zeros(1,time/mesh+1);

% for t = 1:time/mesh+1
%     for s = 1:(pos_bob-pos_alice)/mesh+1
%         % fun(t,s) = j10(abs(psiMinusMinus(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiMinusPlus(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiPlusMinusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2,abs(psiPlusPlusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2)+j01(abs(psiMinusMinus(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiMinusPlus(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2,abs(psiPlusMinusNear(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2,abs(psiPlusPlusNear(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2);
% 
%         % PLOTTING ALICE USING FIRST TWO TERMS OF (31) FROM ARRIVAL TIMES PAPER
%         fun_photon_alice(t,s) = abs(psiPlusMinusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2+abs(psiPlusPlusNear(times(t),pos_alice,times(t),space(s),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang)).^2;
% 
%         % PLOTTING BOB USING SECOND TWO TERMS OF (31) FROM ARRIVAL TIMES PAPER
%         fun_electron_alice(t,s) = abs(psiMinusPlus(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2+abs(psiPlusPlusFar(times(t),space(s),times(t),pos_alice,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega)).^2;
%     end
%     mu_fun_photon_alice(t) = trapz(fun_photon_alice(t,:))*mesh;
%     mu_fun_electron_alice(t) = trapz(fun_electron_alice(t,:))*mesh;
% end


%PLOTTING STATISTICS BELOW
txt = {strcat('theta\_ph=',string(theta_ph)),strcat('theta\_el=',string(theta_el)),strcat('mu\_ph=',string(mu_ph)),strcat('mu\_el=',string(mu_el)),strcat('sigma\_ph=',string(sigma_ph)),strcat('sigma\_el=',string(sigma_el)),strcat('k\_ph=',string(k_ph)),strcat('k\_el=',string(k_el)),strcat('omega=',string(omega))};

for x=1:N
    c = cmap(randi(size(cmap,1)),:);
    figure(1)
    plot(traj_noninteracting_alice(x,:),tt_noninteracting,traj_noninteracting_bob(x,:),tt_noninteracting,'Color',c);
    drawnow
    hold on;

    figure(2)
    plot(traj_interacting_alice(x,:),tt_interacting,traj_interacting_bob(x,:),tt_interacting,'Color',c);
    drawnow
    hold on;
end
hold off;

figure(1)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for noninteracting configuration');
xlim([mu_ph-1,mu_el+1])
ylim([0 time])
text(pos_bob-.5,time-.5,txt);

figure(2)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for interacting configuration');
xlim([mu_ph-1,mu_el+1])
ylim([0 time])
text(pos_bob-.5,time-.5,txt);

figure(3)
scatter(initvals(1,:),yy_photon_alice);
hold on;
scatter(initvals(2,:),yy_electron_alice);
hold on;
xlabel('Initial position','FontSize',20);
ylabel('Arrival time','FontSize',20);
title('Arrival time of photon vs electron at Alice');
xlim([mu_ph-1,mu_el+1])
ylim([0 time])
legend('interacting photon','noninteracting electron');
text(pos_bob-.5,time-.5,txt);

figure(4)
cdfplot(yy_photon_alice);
hold on
cdfplot(yy_electron_alice);
title('Cumulative distribution of arrival times at Alice');
xlim([0 time])
legend('interacting photon','noninteracting electron');
text(time-1,1,txt);

figure(5)
histogram(yy_photon_alice,ceil(N/5));
hold on
histogram(yy_electron_alice,ceil(N/5));
title('Histogram of arrival times at Alice');
xlim([0 time])
legend('interacting photon','noninteracting electron');
text(1,5,txt);

% figure(6)
% xlim([0 time])
% plot(times,mu_fun_photon_alice);
% hold on
% plot(times,mu_fun_electron_alice);
% title('Probability density function of arrival times at Alice')
% legend('interacting photon','noninteracting electron');
% text(1,10,txt);
% 
% figure(7)
% xlim([0 time])
% plot(times,cumtrapz(times,mu_fun_photon_alice));
% hold on
% plot(times,cumtrapz(times,mu_fun_electron_alice));
% title('Cumulative distribution function of arrival times at Alice')
% legend('interacting photon','noninteracting electron');
% text(1,10,txt);

save('interacting_photon_electron_variables.mat')