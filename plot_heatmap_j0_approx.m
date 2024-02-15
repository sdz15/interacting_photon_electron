mu = 0;
sigma = .005;
kappa = -2;
omega = 5;
ang = 0;
time = 2;
mesh = 1/100;

space = (mu-1:mesh:mu+1);
times = (0:mesh:time);

M = struct('cdata',[],'colormap',[]);

for t = 1:time/mesh+1
    frame = NaN(size(space));
    for i = 1:2/mesh+1
        pM = abs(phiMinusApprox(times(t),space(i),sigma,kappa,mu,omega)).^2;
        pP = abs(phiPlusApprox(times(t),space(i),sigma,kappa,mu,omega)).^2;
        frame(i) = j0(pM,pP);
    end
    plot(space,frame)
    xlim([-1 1])
    ylim([0 5])
    xlabel('position')
    M(t) = getframe;
end

% v = VideoWriter('j_00.avi','Motion JPEG AVI');
% open(v);
% writeVideo(v,M);
% close(v);