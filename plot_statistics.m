txt = {strcat('theta\_ph=',string(theta_ph)),strcat('theta\_el=',string(theta_el)),strcat('mu\_ph=',string(mu_ph)),strcat('mu\_el=',string(mu_el)),strcat('sigma\_ph=',string(sigma_ph)),strcat('sigma\_el=',string(sigma_el)),strcat('k\_ph=',string(k_ph)),strcat('k\_el=',string(k_el)),strcat('omega=',string(omega))};

for x=1:N
    figure(1)
    plot(traj_interacting_alice(x,:),tt_interacting,'Color',[0 0.4470 0.7410]);
    drawnow
    hold on;
    plot(traj_single_free(x,:),tt_single_free,'Color',[0.8500 0.3250 0.0980]);
    drawnow
    hold on;
    plot(traj_single_boundary(x,:),tt_single_boundary,'Color',[0.9290 0.6940 0.1250]);
    drawnow
    hold on;
end
hold off;

figure(1)
xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories');
xlim([mu_ph-1,mu_el+1])
ylim([0 time])
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);

figure(2)
scatter(initvals(1,:),yy_photon_alice);
hold on;
scatter(initvals(2,:),yy_electron_alice);
hold on;
scatter(initvals(2,:),yy_single);
hold on;
xlabel('Initial position','FontSize',20);
ylabel('Arrival time','FontSize',20);
title('Arrival time of photon vs electron at Alice');
xlim([mu_ph-1,mu_el+1])
ylim([0 time])
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);

figure(3)
cdfplot(yy_photon);
hold on
cdfplot(yy_single_free);
hold on
cdfplot(yy_single_boundary);
xlabel('Time','Fontsize',20)
title('Cumulative distribution of arrival times at Alice');
xlim([0 time])
ylim([0 1])
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(4)
histogram(yy_photon);
hold on
histogram(yy_single_free);
hold on
histogram(yy_single_boundary);
xlabel('Time','Fontsize',20)
title('Histogram of arrival times at Alice');
xlim([0 time])
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(5)
scatter(initvals(1,:),n_yy_photon);
hold on;
scatter(initvals(2,:),n_yy_single_free);
hold on;
scatter(initvals(2,:),n_yy_single_boundary);
hold on;
xlabel('Initial position','FontSize',20);
ylabel('Arrival time','FontSize',20);
title('Normalized arrival time of photon vs electron at Alice');
xlim([mu_ph-1,mu_el+1])
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)+ylimits(1))*1/2,txt);
hold off

figure(6)
cdfplot(n_yy_photon);
hold on
cdfplot(n_yy_single_free);
hold on
cdfplot(n_yy_single_boundary);
xlabel('Time','Fontsize',20)
title('Normalized cumulative distribution of arrival times at Alice');
xlim([0 time])
ylim([0 1])
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))*3/4,(ylimits(2)+ylimits(1))*1/2,txt);
hold off

figure(7)
histogram(n_yy_photon);
hold on
histogram(n_yy_single_free);
hold on
histogram(n_yy_single_boundary);
xlabel('Time','Fontsize',20)
title('Normalized histogram of arrival times at Alice');
xlim([0 time])
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))*3/4,(ylimits(2)+ylimits(1))*1/2,txt);
hold off

figure(8)
xlim([0 time])
plot(times,mu_fun_photon);
hold on
plot(times,mu_fun_single_free);
hold on
plot(times,mu_fun_single_boundary);
xlabel('Time','Fontsize',20)
title('Probability density function of arrival times at Alice')
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(9)
xlim([0 time])
plot(times,cumtrapz(times,mu_fun_photon));
hold on
plot(times,cumtrapz(times,mu_fun_single_free));
hold on
plot(times,cumtrapz(times,mu_fun_single_boundary));
xlabel('Time','Fontsize',20)
title('Cumulative distribution function of arrival times at Alice')
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off