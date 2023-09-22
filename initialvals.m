% function arr = initialvals(theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el,N,res,mesh,bound)
%     Xin = -bound:mesh:bound;
%     Yin = -bound:mesh:bound;
% 
%     Xmat = ones(length(Yin),1)*Xin;
%     Ymat = Yin'*ones(1,length(Xin));
%     Dist = zeros(2*bound*(1/mesh)+1,2*bound*(1/mesh)+1);
% 
%     for i=1:(2*bound*(1/mesh)+1)
%         for j=1:(2*bound*(1/mesh)+1)
%             psiMMR = abs(psiMinusMinusR((i-(bound*(1/mesh)+1))*mesh,(j-(bound*(1/mesh)+1))*mesh,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el)).^2;
%             psiMPR = abs(psiMinusPlusR((i-(bound*(1/mesh)+1))*mesh,(j-(bound*(1/mesh)+1))*mesh,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el)).^2;
%             psiPMR = abs(psiPlusMinusR((i-(bound*(1/mesh)+1))*mesh,(j-(bound*(1/mesh)+1))*mesh,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el)).^2;
%             psiPPR = abs(psiPlusPlusR((i-(bound*(1/mesh)+1))*mesh,(j-(bound*(1/mesh)+1))*mesh,theta_ph,theta_el,sigma_ph,sigma_el,k_ph,k_el,mu_ph,mu_el)).^2;
%             Dist(i,j) = j00(psiMMR,psiMPR,psiPMR,psiPPR);
%         end
%     end
% 
%     vals=zeros(2,N);
%     for i=1:N
%         [vals(1,i),vals(2,i)]=pinky(Xin,Yin,Dist,res);
%     end
% 
%     arr = vals;
% 
% end

function vals = initialvals(sigma,mu,N,num_particles)
    arr = [];
    for i=1:num_particles
        arr = [arr; repmat(mu(i),1,N)+randn(1,N)*sigma(i)];
    end
    
    vals = arr;

end