% Likelihood function for male translocated inividuals - note this is
% equivalent to a CJS model and we do not estimate N as we are conditioning
% on the first capture (release) of these individuals.

function L = MT(theta,x)

% specify time intervals and time-varying covariates for LNG
time_int = [0.07 0.29 0.05 0.12 0.29 0.25 0.67 0.45 0.91 0.74 1.03 0.39 0.55 0.43 0.55 0.41 1.12 1.39 0.36];
temp = [0.85290000 0.65833333 0.71666667 0.35000000 0.30000000 0.85833333 0.28333333 0.66666667 0.06666667 0.40000000 0.89166667 1.00000000 0.00000000 0.91666667 0.13333333 0.85000000 0.25000000 0.08333333 0.83333333 0.06666667];
% temp = [27.64	26.40	27.28	24.80	24.87	28.10	24.24	27.34	23.71	24.99	28.18	28.12	22.77	28.48	23.95	28.24	24.69	23.81	28.00   22.77];
% temp = rescale(temp);
% moon = [0.75	0.63	0.25	0.75	0.25	0.38	0.25	0.25	0.25	0.50	0.75	0.50	0.50	0.25	1.00	0.88	0.75	0.75    0.00    1.00];
moon = [0.7437361 0.6400749 0.2486548 0.7586744 0.3567393 0.3626744 0.2334367 0.3595274 0.330166 0.4258986 0.8474750 0.5835009 0.5054162 0.2200089 0.8760433 0.9401391 0.6299978 0.7267276 0.1185454 0.8890637];
effort = [2 2 2 2 2 2 2 2 3 3 4 4 4 4 4 4 4 4 4 4];
effort = rescale(effort);

% 

[D,T] = size(x);

% Parameters

% p(t): probability an individual is captured at occasion t given it
% arrived a occasions ago

p(1) = 1;

% % p(.), p(sex), p(origin), p(sex+origin)
% for t = 2:T
%     p(t) = ilogit(theta(33));
% end
 
% % p(time), p(sex+time), p(origin+time), p(origin+sex+time)
% for t = 2:T
%     p(t) = ilogit(theta(31+t));
% end

% p(temp), p(sex+temp), p(origin+temp), p(sex+origin+temp)
for t = 2:T
    p(t) = ilogit(theta(18)+theta(19)*effort(t)+theta(20)*temp(t)+theta(21)*moon(t));
end

% % p(moon), p(sex+moon), p(origin+moon), p(sex+origin+moon)
% for t = 2:T
%     p(t) = ilogit(theta(18)+theta(19)*temp(t));
% end


% phi(t,a): probability an individual in the study at occasion t remians in
% the study until occasion t+1, given that it arrived a occasions ago

% phi(.), phi(sex), phi(origin), phi(origin+sex)
for t = 1:T-1
    phi(t)=(ilogit(theta(22)))^time_int(t);
end

% % phi(time), phi(sex+time), phi(origin+time), phi(origin+sex+time)
% for t = 1:T-1
%     phi(t)=(ilogit(theta(34+t)))^time_int(t);
% end

phi(T)=0;

% find first point and last point
for i = 1:D
    for t = 1:T
        % finds first capture
        % if individual was caught and wasn't caught previously, that time
        % point is added to list
        if x(i,t)==1 && sum(x(i,1:t-1))==0
            f(i)=t;
        end
    end
    % finds last capture
    %for t = T:-1:1
    for t = f(i):T   % correction
        if x(i,t)==1 && sum(x(i,t+1:T))==0
            l(i)=t;
        end
    end
end


% calculate likelihood for each individual and sum
% returns negative likelihood
for i = 1:D
    prob(i) = 0;
    for d = l(i):T
        prodphi1 = 1;
        for j = f(i):d-1
            prodphi1=prodphi1*phi(j);
        end
        prodp1 = 1;
        for j = f(i)+1:d
            prodp1 = prodp1*p(j)^(x(i,j))*(1-p(j))^(1-x(i,j));
        end
        prob(i) = prob(i)+prodphi1*(1-phi(d))*prodp1;
    end
end

prodprob=0;
for i = 1:D
    prodprob = prodprob + log(prob(i));
end

L = -prodprob;