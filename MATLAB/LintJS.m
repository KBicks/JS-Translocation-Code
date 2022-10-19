function Lcombine = LintJS(theta,MWi,FWi)

L1 = JSM(theta,MWi);
L2 = JSF(theta,FWi);

Lcombine = L1+L2;