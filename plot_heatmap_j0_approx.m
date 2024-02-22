mu = 0;
alpha = 1;
kappa = -1;
omega = 5;
ang = 0;
time = 2;
mesh = 1/100;

space = (mu-2:mesh:mu+2);
times = (0:mesh:time);

M = struct('cdata',[],'colormap',[]);

for t = 1:length(times)
    frame = NaN(size(space));
    for i = 1:length(space)
        pM = abs(phiMinusApprox(times(t),space(i),alpha,kappa,mu,omega)).^2;
        pP = abs(phiPlusApprox(times(t),space(i),alpha,kappa,mu,omega)).^2;
        frame(i) = j0(pM,pP);
    end
    plot(space,frame)
    xlim([space(1) space(length(space))])
    ylim([0 2])
    xlabel('position')
    M(t) = getframe;
end