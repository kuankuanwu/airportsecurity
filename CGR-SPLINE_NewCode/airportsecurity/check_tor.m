global c1
global c2
global Beta1
global Beta2
global Budget
global Epsilon
c1 = 2;
c2 = 3;
Beta1 = 5;
Beta2 = 35;
Budget = 95;
Epsilon = 15;
global N
global lambda
global mu1
global mu2
global theta
global RefMatrix
global Lambda
global RefMatrix_tau
RefMatrix_tau = [];
N = 50;
lambda = 40/60;
Lambda = 40/60;
mu1 = 10/60;
mu2 = 1/60;
theta = 0.55;
global L1
global L2
global d1
global d2
d1 = 0.7;
d2 = 0.98;
global chainflag
chainflag=0;
Num_Warmup = 200;
 S1=11;
 S2=18;
%  tau =0.281730242174354;
% tau=0.281578949;
% tau=0.28668511;
global log_test
% iseed =  [271627303,1917625837,3049361867,44444444];
% iseed =  [28532443,984144661,1917997926,44444444];
iseed =  [280402218,686884945,4051340927,44444444];

%  iseed =  [1653381726,2882905302,1955172277,44444444];
 Num_Simulation = 18000 ;
%  Num_Simulation = 18000 ;
Delta0 = 0.001;
log_test_name =['test_check_tor.txt']; 
log_test=fopen(log_test_name, 'w'); 
optimal = 0.362218256;
result=[];


for i = 1 : 10
    RefMatrix_tau = [];
    Tol = 0.1/2^(i-1);
    [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(S1,S2,iseed,...
    Num_Warmup,Num_Simulation,...
    Delta0,Tol);
    result=[result;[Tol,optimal,tau,optimal-tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs]];

end