txt = {strcat('theta\_ph=',string(theta_ph)),strcat('theta\_el=',string(theta_el)),strcat('mu\_ph=',string(mu_ph)),strcat('mu\_el=',string(mu_el)),strcat('sigma\_ph=',string(sigma_ph)),strcat('sigma\_el=',string(sigma_el)),strcat('k\_ph=',string(k_ph)),strcat('k\_el=',string(k_el)),strcat('omega=',string(omega))};

for x=1:N
    figure(1)
    plot(traj_noninteracting_alice(x,:),tt_noninteracting,'Color',[0 0.4470 0.7410]);
    drawnow
    hold on;
    plot(traj_noninteracting_bob(x,:),tt_noninteracting,'Color',[0.8500 0.3250 0.0980]);
    drawnow
    hold on;

    figure(2)
    plot(traj_interacting_alice(x,:),tt_interacting,'Color',[0 0.4470 0.7410]);
    drawnow
    hold on;
    plot(traj_interacting_bob(x,:),tt_interacting,'Color',[0.8500 0.3250 0.0980]);
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
legend({'photon','electron'},'Location','southwest');
text(pos_bob-.5,time-.5,txt);

figure(2)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories for interacting configuration');
xlim([mu_ph-1,mu_el+1])
ylim([0 time])
legend({'photon','electron'},'Location','southwest');
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

figure(6)
xlim([0 time])
plot(times,mu_fun_photon_alice);
hold on
plot(times,mu_fun_electron_alice);
title('Probability density function of arrival times at Alice')
legend('interacting photon','noninteracting electron');
text(1,10,txt);

figure(7)
xlim([0 time])
plot(times,cumtrapz(times,mu_fun_photon_alice));
hold on
plot(times,cumtrapz(times,mu_fun_electron_alice));
title('Cumulative distribution function of arrival times at Alice')
legend('interacting photon','noninteracting electron');
text(1,10,txt);