d1 = 0.7;
d2 = 0.98;
c1 = 3;
c2 = 5;
Beta1 = 5;
Beta2 = 20;
theta = 0.0625;

S1 = 3;
S2 = 4;

ObjFunc = d1*CalcR1(tau,theta) + d2*CalcR2(tau,theta);
WaitingTime = CalcWaitingTime(S1,S2,tau);
Cost = Beta1*CalcP(tau,theta) + Beta2*(1-CalcP(tau,theta)); + S1*c1 + S2*c2;
