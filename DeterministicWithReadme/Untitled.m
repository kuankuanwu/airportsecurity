global c1
global c2
global Beta1
global Beta2
global Budget

c1=1;
c2=2;
Beta1=5;
Beta2=19;
Budget=62;

%(Waiting Time)
global N
global Lambda
global mu1
global mu2
global theta
global Epsilon

Epsilon=30;
N = 50 ;
Lambda=30;  
mu1=4;
mu2=2;
theta=0.25;
ss=[5,19;3,20];
global RefMatrix
RefMatrix=-ones(100,100);
sol=[];
for i=1:2
    [tau, MeanWaitingTime, MinTau_Budget, flag] = GetOptimTau(ss(i,1),ss(i,2),1);
    sol=[sol;tau];
end