% SINGLE ELECTRON

function val = velocity_single_electron(t,s,theta,sigma,k,mu,omega,L,step)
    psiM = abs(psiMinus(t,s,theta,sigma,k,mu,omega,L,step)).^2;
    psiP = abs(psiPlus(t,s,theta,sigma,k,mu,omega,L,step)).^2;
    j_0 = j0(psiM,psiP);
    if (norm(j_0)>=1e-6)
        j_1 = j1(psiM,psiP);
        val = j_1./j_0;
        return
    end
    val = 0;
end