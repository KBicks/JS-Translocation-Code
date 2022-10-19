function Nt_cons = calcNtJS(phi,beta,N,NT)

[nt_start,K] = size(NT);

Nt_cons(1) = beta(1)*N;

for t = 1:K-1
    eta = Nt_cons(t)*phi(t)+N*beta(t+1);
    Nt_cons = [Nt_cons eta];
end
    