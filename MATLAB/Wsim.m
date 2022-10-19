% JS likelihood for wild born male individuals

function L = Wsim(theta,x,time_int,tau)

[D,T] = size(x);
% Parameters

% N: population size (super-population); n denotes the number of uncaught
% individuals and is what we estimate.

n = exp(theta(1));
N = n+D;

% beta(t): proportion of the (N-nt) individuals who are first available for capture at
% occasion t.  Note beta starts at beta(1) here rather than beta(0).  We assume new arrivals can only start from time tau.

for t = 1:tau-1
    beta(t) = 0;
end

sumbeta = 0;
ind = 1;

for t = tau+1:T
    ind = ind+1;
   sumbeta = sumbeta+exp(theta(ind));
end

ind = 1;
for t = tau+1:T
    ind = ind+1;
    beta(t) = exp(theta(ind))/(1+sumbeta);
end

beta(tau) = 1/(1+sumbeta);


% % beta(.)
% beta(tau) = 1/(1+(T-tau)*exp(theta(2)));
% 
% for t = tau+1:T
%     beta(t)= exp(theta(2))/(1+(T-tau)*exp(theta(2)));
% end
    

% p(t): probability an individual is captured at occasion t given it
% arrived a occasions ago

p(1) = 1;

% p(.), p(sex), p(origin), p(sex+origin)
for t = 2:T
    p(t) = ilogit(theta(11));
end
 
% phi(t,a): probability an individual in the study at occasion t remians in
% the study until occasion t+1, given that it arrived a occasions ago

% phi(.), phi(sex), phi(origin), phi(origin+sex)
for t = 1:T-1
    phi(t)=(ilogit(theta(12)))^time_int(t);
end

phi(T)=0;


for i = 1:D
    for t = 1:T
        if x(i,t)==1 && sum(x(i,1:t-1))==0
            f(i)=t;
        end
    end
    for t = T:-1:1
        if x(i,t)==1 && sum(x(i,t+1:T))==0
            l(i)=t;
        end
    end
end




% probabilities for new entrant observed individual i

for i = 1:D
    prob(i) = 0;
    for b = tau:f(i)
        for d = l(i):T
            prodphi1 = 1;
            for j = b:d-1
                prodphi1=prodphi1*phi(j);
            end
            prodp1 = 1;
            for j = b:d
                prodp1 = prodp1*p(j)^(x(i,j))*(1-p(j))^(1-x(i,j));
            end
            prob(i) = prob(i)+beta(b)*prodphi1*(1-phi(d))*prodp1;
        end
    end
end

% probability for unobserved individuals

prob0 = 0;
for b = tau:T
    for d = b:T
        prodphi2 = 1;
        for j = b:d-1
            prodphi2=prodphi2*phi(j);
        end
        prodp2 = 1;
        for j = b:d
            prodp2 = prodp2*(1-p(j));
        end
        prob0 = prob0+beta(b)*prodphi2*(1-phi(d))*prodp2;
    end
end

prodprob=0;
for i = 1:D
    prodprob = prodprob + log(prob(i));
end

L = gammaln(N+1) - gammaln(n+1) + prodprob + n*log(prob0);

L = -L;