clear all

% set optimisation rules - no need to change these

opt=optimset('Largescale','off', 'TolFun',.000001, 'TolX',.000001, 'MaxIter',100000, 'MaxFunEvals',1000000);

% set capture histories to 4 matrices - male and female translocated and
% wildborn
%[MTr,FTr,MWi,FWi] = LNGdata;

MTr = table2array(readtable("IM_lng_mt.csv"));
FTr = table2array(readtable("IM_lng_ft.csv"));
MWi = table2array(readtable("IM_lng_mw.csv"));
FWi = table2array(readtable("IM_lng_fw.csv"));

% fit model using separate likelihood function for each group and summing
% likelihoods
[y,yval,exitflag,output,grad,hessian] = fminunc(@Lint,[4,4,rand(1,20)],opt,MTr,FTr,MWi,FWi);
[z,zval,exitflagJS,outputJS,gradJS,hessianJS] = fminunc(@LintJS,[4,4,rand(1,40)],opt,[MTr;MWi],[FTr;FWi]);

[NM,NF,betaM,betaF,phiM,phiF,pM,pF]  = paramcalc(y,MWi,FWi);
[NMJS,NFJS,betaMJS,betaFJS,phiMJS,phiFJS,pMJS,pFJS]  = paramcalcJS(z,[MTr;MWi],[FTr;FWi]);

% Number of females over time
% NtF = calcNt(phiF,betaF,NF,FTr);
NtFJS = calcNtJS(phiFJS,betaFJS,NFJS,FTr);

% Number of males over time
% NtM = calcNt(phiM,betaM,NM,MTr);
NtMJS = calcNtJS(phiMJS,betaMJS,NMJS,MTr);

% Bootstrap for new model - note the number of iterations will need to be
% increased to ~250

for s = 1:250
    BMTr = datasample(MTr,size(MTr,1),1);
    BFTr = datasample(FTr,size(FTr,1),1);
    BMWi = datasample(MWi,size(MWi,1),1);
    BFWi = datasample(FWi,size(FWi,1),1);
    
%     y(s,:) = fminunc(@Lint,[4,4,rand(1,20)],opt,BMTr,BFTr,BMWi,BFWi);
%     [NMb(s),NFb(s),betaMb(s,:),betaFb(s,:),phiMb(s,:),phiFb(s,:),pMTb(s,:),pFTb(s,:),pMWb(s,:),pFWb(s,:)]  = paramcalc(y(s,:),BMWi,BFWi);

    z(s,:) = fminunc(@LintJS,[4,4,rand(1,40)],opt,[BMTr;BMWi],[BFTr;BFWi]);
    [NMJSb(s),NFJSb(s),betaMJSb(s,:),betaFJSb(s,:),phiMJSb(s,:),phiFJSb(s,:),pMJSb(s,:),pFJSb(s,:)]  = paramcalcJS(z(s,:),[BMTr;BMWi],[BFTr;BFWi]);

% Number of females over time
% NtFb(s,:) = calcNt(phiFb(s,:),betaFb(s,:),NFb(s),BFTr);
NtFJSb(s,:) = calcNtJS(phiFJSb(s,:),betaFJSb(s,:),NFJSb(s),BFTr);

% Number of males over time
% NtMb(s,:) = calcNt(phiMb(s,:),betaMb(s,:),NMb(s),BMTr);
NtMJSb(s,:) = calcNtJS(phiMJSb(s,:),betaMJSb(s,:),NMJSb(s),BMTr);



end


% % Fitting the standard JS model to the LNG data
% z = fminunc(@LintJS,[4,4,zeros(1,52)],opt,MTr,FTr,MWi,FWi);
% 
% [NMJS,NFJS,betaMJS,betaFJS,phiMJS,phiFJS,pMJS,pFJS]  = paramcalcJS(z,MTr,FTr,MWi,FWi);
% 
% % Number of females over time
% NtFJS = calcNtJS(phiFJS,betaFJS,NFJS,FTr);
% 
% % Number of males over time
% NtMJS = calcNtJS(phiMJS,betaMJS,NMJS,MTr);

% Bootstrap for JS model

% for s = 1:3
%     BMTr = datasample(MTr,size(MTr,1),1);
%     BFTr = datasample(FTr,size(FTr,1),1);
%     BMWi = datasample(MWi,size(MWi,1),1);
%     BFWi = datasample(FWi,size(FWi,1),1);
%     
%     z(s,:) = fminunc(@LintJS,[4,4,rand(1,52)],opt,BMTr,BFTr,BMWi,BFWi);
%     [NMJSb(s),NFJSb(s),betaMJSb(s,:),betaFJSb(s,:),phiMJSb(s,:),phiFJSb(s,:),pMJSb(s,:),pFJSb(s,:)]  = paramcalcJS(z(s,:),BMTr,BFTr,BMWi,BFWi);
% 
% % Number of females over time
% NtFJSb(s,:) = calcNtJS(phiFJSb(s,:),betaFJSb(s,:),NFJSb(s),BFTr);
% 
% % Number of males over time
% NtMJSb(s,:) = calcNtJS(phiMJSb(s,:),betaMJSb(s,:),NMJSb(s),BMTr);
% end


