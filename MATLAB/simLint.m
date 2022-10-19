function Lcombine = simLint(theta,Tr,Wi,time_intervals,tau)

L1 = Tsim(theta,Tr,time_intervals);
L2 = Wsim(theta,Wi,time_intervals,tau);

Lcombine = L1+L2;