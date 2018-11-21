global N
global Lambda
global mu1
global mu2
global theta
global Epsilon
N = 50;
Lambda = 2;
mu1 = 0.66667;
mu2 = 0.016667;
theta = 0.55;
Epsilon = 25;

RefMatrix = ones(100,100);
mm=zeros(30,1);
tau=linspace(0,1);
% [AA] = xlsread('aa','1','A1:C30');
% for j = 1:30
%[WaitingTime, PE1, PE2, PE12] = MeanWait(AA(i,3),AA(i,1),AA(i,2));
%[tau,MeanWaitingTime,MinTau_Budget,flag] = GetOptimTau(AA(j,1),AA(j,2),1);
[WaitingTime, PE1, PE2, PE12] = MeanWait_MatrixInverse(0.1,3,20);
% mm(j)=WaitingTime*60-Epsilon;
% end
% xlswrite('mm',mm);
plot(tau,WaitingTime);