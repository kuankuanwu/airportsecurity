
% InputFileName = 'input1.txt';
% InitializeAndReadInput(InputFileName)

%(Constraint)
% addpath('airportsecurity');  
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
%(Waiting Time)
global N
global Lambda
global mu1
global mu2
global theta
N=50;
Lambda = 40/60;
mu1 = 10/60;
mu2 = 1/60;
theta = 0.55;
global L1
global L2

global chainflag
chainflag=0;
% S1=5;
% S2=5;
% tau=0.200363790988488;
res=[];
  
global W1
global W2
% iseed      = [11111111,22222222,33333333,44444444];  
% Num_Warmup=200;
% Num_Simulation=81000;

S1=2;
S2=23;
tau=0.24803;
P = CalcP(tau, theta); 
RefMatrix=-ones(100,100);    
[MeanWaitingTime, PE1, PE2, PE12] = MeanWait_MatrixInverse(tau,S1,S2);

% for i =1:30
% S1=x_best(i,3);
% S2=x_best(i,4);
% tau=x_best(i,5);
% P = CalcP(tau, theta); 
% RefMatrix=-ones(100,100);    
% [MeanWaitingTime, PE1, PE2, PE12] = MeanWait_MatrixInverse(tau,S1,S2);
% % MeanWaitingTime   = MeanWait_Chain(tau,S1,S2);
% % [tau, MeanWaitingTime, MinTau_Budget, flag] = GetOptimTau(S1,S2,1);
% 
% res=[res;[i,MeanWaitingTime-Epsilon]];
% % disp(['MeanWaitingTime=',num2str(MeanWaitingTime)])
% % res=[res;[Lambda,P,L1,L2,PE1,PE2,PE12,WaitingTime]];
% %  res=[res;[Lambda,L1,L2,W1,W2,MeanWaitingTime,tau]];
% end
% S=[23,45;34,34];
% ans_out=[];
% for i=1:2
%     S1=S(i,1);S2=S(i,2); 
%     MinTau_Budget = GetMinTau(S1,S2);
%     [tau , MeanWaitingTime, flag] = SolveTau(Epsilon,S1,S2,MinTau_Budget);
%     ans_out=[ans_out;[S1,S2,tau,MeanWaitingTime]];
% end