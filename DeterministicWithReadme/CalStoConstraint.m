% InputFileName = ['input',num2str(1),'.txt'];
% InitializeAndReadInput(InputFileName)
addpath('D:\Google Drive\碩士\project\航空站問題\CGR-SPLINE_NewCode\airportsecurity')

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
Num_Warmup = 200;
% iseed      = [1589895398,166326175,2608342239,44444444];
% 
h1=[];
h2=[];
res=[];

%    Num_Simulation = 81000 ;
%    S1=8;
%    S2=4;
%    tau = 0.252220342478897;
%    [MeanWaitingTime, PE1, PE2, PE12]=MeanWait_MatrixInverse(tau,S1,S2);
%    max_securitylevel = d1*CalcR1(tau,theta) + d2*CalcR2(tau,theta);
%    Cost = c1*S1 + c2*S2 + Beta1*CalcP(tau,theta) + Beta2*(1-CalcP(tau,theta));
%    [MeanWaitingTimeSim,SE]=Simulation_AirportModel(tau,S1,S2,iseed,Num_Warmup,Num_Simulation);
% 
%    [MeanWaitingTime, PE1, PE2, PE12]=MeanWait_MatrixInverse(tau,S1,S2);
%   h1=[h1;[S1,S2,tau,max_securitylevel,MeanWaitingTimeSim-Epsilon,Cost-Budget,MeanWaitingTime-Epsilon]];

for i=1:25
%    test = [];
%    X=Abest(i).x;
%    S1=4;
%    S2=22;
%    tau=0.281730242;
%    iseed =  [1653381726,2882905302,1955172277,44444444];
%    Num_Simulation = 54000 ;
%    [MeanWaitingTime, PE1, PE2, PE12]=MeanWait_MatrixInverse(tau,S1,S2);
%    max_securitylevel = d1*CalcR1(tau,theta) + d2*CalcR2(tau,theta);
%    Cost = c1*S1 + c2*S2 + Beta1*CalcP(tau,theta) + Beta2*(1-CalcP(tau,theta));
%    [MeanWaitingTimeSim,SE]=Simulation_AirportModel(tau,S1,S2,iseed,Num_Warmup,Num_Simulation);
%    test=[test;[S1,S2,tau,max_securitylevel,MeanWaitingTimeSim-Epsilon,Cost-Budget,SE,Awork(i),MeanWaitingTime-Epsilon,Acurr(i).lastRAipseed(1),Acurr(i).lastRAipseed(2),Acurr(i).lastRAipseed(3),Acurr(i).lastRAipseed(4)]];

   X=Acurr(i).x;
   S1=Acurr(i).x(1);
   S2=Acurr(i).x(2);
   tau=Acurr(i).fn;
   iseed = Acurr(i).lastRAipseed;
   Num_Simulation = Acurr(i).mk ;
   [MeanWaitingTime, PE1, PE2, PE12]=MeanWait_MatrixInverse(tau,S1,S2);
   max_securitylevel = d1*CalcR1(tau,theta) + d2*CalcR2(tau,theta);
   Cost = c1*S1 + c2*S2 + Beta1*CalcP(tau,theta) + Beta2*(1-CalcP(tau,theta));
   [MeanWaitingTimeSim,SE]=Simulation_AirportModel(tau,S1,S2,iseed,Num_Warmup,Num_Simulation);
   h1=[h1;[S1,S2,tau,max_securitylevel,MeanWaitingTimeSim-Epsilon,Cost-Budget,SE,Awork(i),MeanWaitingTime-Epsilon,Acurr(i).lastRAipseed(1),Acurr(i).lastRAipseed(2),Acurr(i).lastRAipseed(3),Acurr(i).lastRAipseed(4)]];

   Num_Simulation = Abest(i).mk ;
   S1=Abest(i).x(1);
   S2=Abest(i).x(2);
   iseed = Abest(i).lastRAipseed;
   tau = Abest(i).fn;
    [MeanWaitingTime, PE1, PE2, PE12]=MeanWait_MatrixInverse(tau,S1,S2);
    max_securitylevel = d1*CalcR1(tau,theta) + d2*CalcR2(tau,theta);
    Cost = c1*S1 + c2*S2 + Beta1*CalcP(tau,theta) + Beta2*(1-CalcP(tau,theta));
    [MeanWaitingTimeSim,SE]=Simulation_AirportModel(tau,S1,S2,iseed,Num_Warmup,Num_Simulation);
    h2=[h2 ;[S1,S2,tau,max_securitylevel,MeanWaitingTimeSim-Epsilon,Cost-Budget,SE,Awork(i),MeanWaitingTime-Epsilon,Abest(i).mk,Abest(i).lastRAipseed(1),Abest(i).lastRAipseed(2),Abest(i).lastRAipseed(3),Abest(i).lastRAipseed(4)]];


end

      