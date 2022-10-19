% Likelihood function for male translocated inividuals - note this is
% equivalent to a CJS model and we do not estimate N as we are conditioning
% on the first capture (release) of these individuals.

function L = Tsim(theta,x,time_int)

[D,T] = size(x);

% Parameters

% p(t): probability an individual is captured at occasion t given it
% arrived a occasions ago

p(1) = 1;

% p(.), p(sex), p(origin), p(sex+origin)
for t = 2:T
    p(t) = ilogit(theta(11));
end
 
% phi(.), phi(sex), phi(origin), phi(origin+sex)
for t = 1:T-1
    phi(t)=(ilogit(theta(12)))^time_int(t);
end

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