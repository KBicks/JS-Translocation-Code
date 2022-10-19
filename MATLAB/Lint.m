function Lcombine = Lint(theta,MTr,FTr,MWi,FWi)

L1 = MT(theta,MTr);
L2 = FT(theta,FTr);
L3 = MW(theta,MWi);
L4 = FW(theta,FWi);

Lcombine = L1+L2+L3+L4;