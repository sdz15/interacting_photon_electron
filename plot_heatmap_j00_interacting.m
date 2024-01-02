theta_ph = 0;
theta_el = pi/2;
mu_ph = 0;
mu_el = 1;
sigma_ph = .005;
sigma_el = .005;
k_ph = 0;
k_el = 0;
omega = 2;
ang = 0;
time = 4;
N = 500;
mesh = .01;
step = 10/mesh;
pos_alice = mu_ph-1;
pos_bob = mu_el+1;
L = pos_alice-mu_el;

space = (mu_ph-1:mesh:mu_el+1);
times = (0:mesh:time);

M = struct('cdata',[],'colormap',[]);

for t = 1:time/mesh+1
    frame = NaN(size(space));
    for i = 1:(mu_el-mu_ph+2)/mesh+1
        for j = 1:(mu_el-mu_ph+2)/mesh+1
            if (space(j)>=space(i))
                pMM = abs(psiMinusMinus(times(t),space(i),times(t),space(j),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
                pMP = abs(psiMinusPlus(times(t),space(i),times(t),space(j),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
                if (space(i)+times(t)-1e-10 > space(j)-times(t))
                    pPM = abs(psiPlusMinusNear(times(t),space(i),times(t),space(j),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;
                    pPP = abs(psiPlusPlusNear(times(t),space(i),times(t),space(j),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,ang,step)).^2;
                else
                    pPM = abs(psiPlusMinusFar(times(t),space(i),times(t),space(j),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
                    pPP = abs(psiPlusPlusFar(times(t),space(i),times(t),space(j),theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,omega,step)).^2;
                end
                frame(i,j) = j00(pMM,pMP,pPM,pPP);
            end
        end
    end
    surf(space,space,frame);
    zlim([0 1])
    xlabel('photon')
    ylabel('electron')
    zlabel('j00')
    M(t) = getframe;
end
movie(M,100)

% v = VideoWriter('j_00.avi','Motion JPEG AVI');
% open(v);
% writeVideo(v,M);
% close(v);