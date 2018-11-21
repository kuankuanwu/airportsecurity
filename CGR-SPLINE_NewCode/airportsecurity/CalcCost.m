global c1
global c2
global Beta1
global Beta2
global Budget
global Epsilon
c1=2;
c2=3;
Beta1=5;
Beta2=35;
Budget=95;
%(Waiting Time)
global N
N=50;
global lambda
global mu1
global mu2
global theta
global cost

Epsilon=15;
lambda=40/60;
mu1=10/60;
mu2=1/60;
theta = 0.55;

S1=6;
S2=21;
tau=0.30196634;

global P
P = CalcP(tau,theta)
cost= Beta1*P+Beta2*(1-P)+c1*S1+c2*S2
