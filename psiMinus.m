function val = psiMinus(t,s,theta,alpha,k,mu,omega,L,step)

if t>=L+s
    f = -phiMinus(t-L-s,-L,theta,alpha,k,mu,omega,step);
else
    f = 0;
end

n = 0;
chisum = 0;

if t > 0
    while (t>(2*n*L))
        chisum = chisum + chi(t,s,theta,alpha,k,mu,omega,L,step,-1,n);
        n = n+1;
    end
end

val = phiMinus(t,s,theta,alpha,k,mu,omega,step)+f+chisum;