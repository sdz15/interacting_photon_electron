theta = pi/2;
mu = 0;
sigma = .005;
k = -5;
omega = 5;
ang = 0;
time = 2;
mesh = 1/100;
step = 1/mesh;

space = (mu-2:mesh:mu+2);
times = (0:mesh:time);

M = struct('cdata',[],'colormap',[]);

for t = 1:time/mesh+1
    frame = NaN(size(space));
    for i = 1:2/mesh+1
        pM = abs(phiMinus(times(t),space(i),theta,sigma,k,mu,omega,step)).^2;
        pP = abs(phiPlus(times(t),space(i),theta,sigma,k,mu,omega,step)).^2;
        frame(i) = j0(pM,pP);
    end
    plot(space,frame)
    xlim([mu-2 mu+2])
    ylim([0 2])
    xlabel('position')
    M(t) = getframe;
end