theta = pi/2;
mu = 0;
alpha = 1;
kappa = 1;
omega = 20;
k = kappa*omega;
ang = 0;
time = 10;
mesh = 1/100;
step = 1/mesh;

space = (-time:mesh:time);
times = (0:mesh:time);

M = struct('cdata',[],'colormap',[]);

for t = 1:length(times)
    frame = NaN(size(space));
    for i = 1:length(space)
        pM = abs(phiMinus(times(t),space(i),theta,alpha,k,mu,omega,step)).^2;
        pP = abs(phiPlus(times(t),space(i),theta,alpha,k,mu,omega,step)).^2;
        pMApprox = abs(phiMinusApprox(times(t),space(i),alpha,kappa,mu,omega)).^2;
        pPApprox = abs(phiPlusApprox(times(t),space(i),alpha,kappa,mu,omega)).^2;
        frame(i) = (j1(pM,pP)/j0(pM,pP))-(j1(pMApprox,pPApprox)/j0(pMApprox,pPApprox));
    end
    plot(space,frame)
    xlim([space(1) space(length(space))])
    ylim([-2 2])
    xlabel('position')
    ylabel('error')
    title({strcat('time=',string(times(t)))})
    M(t) = getframe;
end