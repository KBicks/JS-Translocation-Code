clear all

% set optimisation rules - no need to change these

opt=optimset('Largescale','off', 'TolFun',.000001, 'TolX',.000001, 'MaxIter',100000, 'MaxFunEvals',1000000);

fileloc = 'SIM\';

filesT = dir('SIM\s1_simT*.csv');
filesW = dir('SIM\s1_simW*.csv');

tauT = 4;
time_intervals = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];


for i = 1:250
    
    % load in simulated data files
    T = fullfile(fileloc, filesT(i).name);
    simT = csvread(T,1,0);
    
    W = fullfile(fileloc, filesW(i).name);
    simW = csvread(W,1,0);
    
    % run new model
    y(i,:) = fminunc(@simLint,[2,rand(1,11)],opt,simT,simW,time_intervals, tauT);
    [nT(i), NT(i), betaT(i,:), phiT(i,:), pT(i,:)] = simparamcalc(y(i,:),simW, time_intervals,tauT);

    % run standard JS model
    z(i,:) = fminunc(@simJS,[2,rand(1,11)],opt,[simT;simW],time_intervals);
    [nJS(i), NJS(i), betaJS(i,:), phiJS(i,:), pJS(i,:)] = simparamcalcJS(z(i,:),[simT;simW],time_intervals, 1);
    
    
    
    % Number of males over time
    NtT(i,:) = calcNt(phiT(i,:),betaT(i,:),NT(i),simT);
    NtJS(i,:) = calcNtJS(phiJS(i,:),betaJS(i,:),NJS(i),simT);
    
    
end


    

