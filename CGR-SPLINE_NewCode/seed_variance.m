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

 global log_test
 log_test = "test";

Num_Warmup = 200;
 S1=4;
 S2=22;

 global Num_Simulation
 Num_Simulation = 162000 ;
 tor = 0.25/(Num_Simulation)^0.25;
 tor_vector = 0.001:0.001:0.05;
 
 rng('shuffle'); 
 record = [];
 global RefMatrix_tau
 RefMatrix_tau = [];
for i = 1:length(tor_vector)
    RefMatrix_tau = [];
%     iseed = randi([0,99999999],1,4);
    iseed = [3861239094	3139040435	811603038	44444444];
    tor = tor_vector(i);
    [tau,AvgWaitingTime,Se] = GetOptimTau_Simulation_Bisection(S1,S2,iseed,Num_Warmup,Num_Simulation,0.001,tor) ;
    record = [record; iseed(1),iseed(2),iseed(3),iseed(4),tau,AvgWaitingTime,Se];
end
    

