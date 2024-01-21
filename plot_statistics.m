txt = {strcat('theta\_ph=',string(theta_ph)),strcat('theta\_el=',string(theta_el)),strcat('mu\_ph=',string(mu_ph)),strcat('mu\_el=',string(mu_el)),strcat('sigma\_ph=',string(sigma_ph)),strcat('sigma\_el=',string(sigma_el)),strcat('k\_ph=',string(k_ph)),strcat('k\_el=',string(k_el)),strcat('omega=',string(omega))};

figure(1)
hold on;
xlim([mu_ph-1,pos_bob])
ylim([0 time])


for x=1:N
    figure(1)
    plot(traj_interacting_photon(x,:),tt_interacting_photon,'Color',[0 0.4470 0.7410]);
    drawnow

    plot(traj_single_free(x,:),tt_single_free,'Color',[0.8500 0.3250 0.0980]);
    drawnow
    plot(traj_single_boundary(x,:),tt_single_boundary,'Color',[0.9290 0.6940 0.1250]);
    drawnow
end

xlabel('Position','FontSize',20);
ylabel('Time','FontSize',20);
title('Trajectories');
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*1/4,txt);
hold off

figure(2)
hold on;
scatter(initvals(1,:),yy_photon);
scatter(initvals(2,:),yy_single_free);
scatter(initvals(2,:),yy_single_boundary);
xlabel('Initial position','FontSize',20);
ylabel('Arrival time','FontSize',20);
title('Arrival time of photon vs electron at Alice');
xlim([mu_ph-1,mu_el+1])
ylim([0 time])
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(3)
hold on;
cdfplot(yy_photon);
cdfplot(yy_single_free);
cdfplot(yy_single_boundary);
% plot(times_prob,cumtrapz(times_prob,mu_fun_photon));
plot(times_prob,cumtrapz(times_prob,mu_fun_single_free));
plot(times_prob,cumtrapz(times_prob,mu_fun_single_boundary));
xlabel('Time','Fontsize',20)
title('Cumulative distributions of arrival times at Alice');
xlim([0 time])
ylim([0 1])
legend({'interacting photon','noninteracting electron','single electron w/ boundary','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(4)
hold on;
histogram(yy_photon);
histogram(yy_single_free);
histogram(yy_single_boundary);
% plot(times_prob,mu_fun_photon);
plot(times_prob,10*mu_fun_single_free);
plot(times_prob,10*mu_fun_single_boundary);
xlabel('Time','Fontsize',20)
title('Histogram and pdf of arrival times at Alice');
xlim([0 time])
legend({'interacting photon','noninteracting electron','single electron w/ boundary','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off

figure(5)
hold on;
% scatter(initvals(1,:),n_yy_photon);
% scatter(initvals(2,:),n_yy_single_free);
% scatter(initvals(2,:),n_yy_single_boundary);
% xlabel('Initial position','FontSize',20);
% ylabel('Arrival time','FontSize',20);
% title('Normalized arrival time of photon vs electron at Alice');
% xlim([mu_ph-1,mu_el+1])
% legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
% xlimits=xlim;
% ylimits=ylim;
% text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)+ylimits(1))*1/2,txt);
scatter(initvals(1,:),yy_photon);
scatter(initvals(2,:)-mu_el,yy_single_free);
scatter(initvals(2,:)-mu_el,yy_single_boundary);
xlabel('Initial position','FontSize',20);
ylabel('Arrival time','FontSize',20);
title('Arrival times of photon vs electron at Alice, Bob shifted to mu\_ph');
xlim([mu_ph-1,mu_el+1])
ylim([0 time])
legend({'interacting photon','noninteracting electron','single electron w/ boundary'},'Location','southwest');
xlimits=xlim;
ylimits=ylim;
text(xlimits(1)+(xlimits(2)-xlimits(1))/16,(ylimits(2)-ylimits(1))*5/8,txt);
hold off