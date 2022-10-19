function Nt_cons = calcNt(phi,beta,N,DT)

[nT,T]=size(DT);
Nt_cons(1)=nT;

for t = 1:T-1
    eta = Nt_cons(t)*phi(t)+N*beta(t+1);
    Nt_cons = [Nt_cons eta];
end
    