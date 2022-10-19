function [NM,NF,betaM,betaF,phiM,phiF,pMT,pFT,pMW,pFW] = paramcalc(theta,MWi,FWi)

% specify time intervals and time-varying covariates for LNG
time_int = [0.07 0.29 0.05 0.12 0.29 0.25 0.67 0.45 0.91 0.74 1.03 0.39 0.55 0.43 0.55 0.41 1.12 1.39 0.36];
temp = [0.85290000 0.65833333 0.71666667 0.35000000 0.30000000 0.85833333 0.28333333 0.66666667 0.06666667 0.40000000 0.89166667 1.00000000 0.00000000 0.91666667 0.13333333 0.85000000 0.25000000 0.08333333 0.83333333 0.06666667];
% temp = [27.64	26.40	27.28	24.80	24.87	28.10	24.24	27.34	23.71	24.99	28.18	28.12	22.77	28.48	23.95	28.24	24.69	23.81	28.00   22.77];
% temp = rescale(temp);
% moon = [0.75	0.63	0.25	0.75	0.25	0.38	0.25	0.25	0.25	0.50	0.75	0.50	0.50	0.25	1.00	0.88	0.75	0.75    0.00    1.00];
moon = [0.7437361 0.6400749 0.2486548 0.7586744 0.3567393 0.3626744 0.2334367 0.3595274 0.330166 0.4258986 0.8474750 0.5835009 0.5054162 0.2200089 0.8760433 0.9401391 0.6299978 0.7267276 0.1185454 0.8890637];
effort = [2 2 2 2 2 2 2 2 3 3 4 4 4 4 4 4 4 4 4 4];
effort = rescale(effort);

Dm = size(MWi,1);
Df = size(FWi,1);
T = size(MWi,2);

tau = 5;

% Parameters

% Nm: male population size (super-population);

nM = exp(theta(1));
NM = nM+Dm;

% Nf: female population size (super-population);

nF = exp(theta(2));
NF = nF+Df;

% betaM(t)/betaF(t): proportion of the (N-nt) individuals who are first available for capture at
% occasion t.  Note beta starts at beta(1) here rather than beta(0).  We assume new arrivals can only start from time tau.

for t = 1:tau-1
    betaM(t) = 0;
end

sumbeta = 0;
ind = 2;

for t = tau+1:T
    ind = ind+1;
   sumbeta = sumbeta+exp(theta(ind));
end

ind = 2;
for t = tau+1:T
    ind = ind+1;
    betaM(t) = exp(theta(ind))/(1+sumbeta);
end

betaM(tau) = 1/(1+sumbeta);

% for beta(t) - below for beta(t*sex)
betaF = betaM;

% for t = 1:tau-1
%     betaF(t) = 0;
% end
% 
% sumbeta = 0;
% ind = 16;
% 
% for t = tau+1:T
%     ind = ind+1;
%    sumbeta = sumbeta+exp(theta(ind));
% end
% 
% ind = 16;
% for t = tau+1:T
%     ind = ind+1;
%     betaF(t) = exp(theta(ind))/(1+sumbeta);
% end
% 
% betaF(tau) = 1/(1+sumbeta);

% pM(t)/pF(t): probability a male/female individual is captured at occasion t given it
% arrived a occasions ago
    
pMT(1) = 1;
% p(temp)male
for t = 2:T
    pMT(t) = ilogit(theta(18)+theta(19)*effort(t)+theta(20)*temp(t)+theta(21)*moon(t));
end

pFT = pMT;
pMW = pMT;
%
% pF(1) = 1;
% % p(temp) female
% for t = 2:T
%     pF(t) = ilogit(theta(20)+theta(21)*temp(t));
% end

% for t = 1:tau-1
%     pMW(t) = 0;
% end
% % p(temp)male
% for t = tau:T
%     pMW(t) = ilogit(theta(22)+theta(23)*effort(t)+theta(24)*temp(t)+theta(25)*moon(t));
% end

pFW = pMW;


% phiM(t)/F(t): probability an individual in the study at occasion t remians in
% the study until occasion t+1

for t = 1:T-1
    phiM(t)=(ilogit(theta(22)))^time_int(t);
end

for t = 1:T-1
    phiF(t)=(ilogit(theta(22)))^time_int(t);
end


