function val = R(j,s,xi,omega,L,n)
x1 = -omega*besselfun(radical(xi,(2*n+1)*L,s),omega,n+j+1).*((2*n+1)*L-s).*((xi-((2*n+1)*L-s)).^(n+j));
x2 = 0;

if (n+j>0)
    x2 = (n+j).*besselfun(radical(xi,(2*n+1)*L,s),omega,n+j).*((xi-((2*n+1)*L-s)).^(n+j-1));
end

val = x1+x2;