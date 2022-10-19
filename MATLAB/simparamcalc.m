function [n,N,beta,phi,p] = simparamcalc(theta,MWi,time_int,tau)

Dm = size(MWi,1);
T = size(MWi,2);

% Parameters

% Nm: male population size (super-population);

n = exp(theta(1));
N = n+Dm;


% betaM(t)/betaF(t): proportion of the (N-nt) individuals who are first available for capture at
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



% beta(tau) = 1/(1+(T-tau)*exp(theta(2)));
% 
% for t = tau+1:T
%     beta(t)= exp(theta(2))/(1+(T-tau)*exp(theta(2)));
% end


% pM(t)/pF(t): probability a male/female individual is captured at occasion t given it
% arrived a occasions ago
    
p(1) = 1;
% p(temp)male
for t = 2:T
    p(t) = ilogit(theta(11));
end

% phiM(t)/F(t): probability an individual in the study at occasion t remians in
% the study until occasion t+1

for t = 1:T-1
    phi(t)=(ilogit(theta(12)))^time_int(t);
end



