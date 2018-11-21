%  global N
%     global Lambda
%     global mu1
%     global mu2
%     global theta
% global Epsilon
% N = 50;
% Lambda = 20;
% mu1 = 3;
% mu2 = 2;
% theta = 0.25;

% tau=linspace(0,1);
% meantime=MeanWait_MatrixInverse(tau,3,20);
% plot(tau,meantime);

clc
global chainflag
chainflag=0;
global c1
global c2
global Beta1
global Beta2
global Budget
global Epsilon
c1=3;
c2=5;
Beta1=5;
Beta2=20;
Budget=44;
%(Waiting Time)
global N
N=50;
global Lambda
global mu1
global mu2
global theta


Epsilon=1;
Lambda=18;
mu1=8;
mu2=6;
theta = 0.0625;
global L1
global L2

% InputFileName = 'input1.txt';
% InitializeAndReadInput(InputFileName)


S1=4; 
S2=4;

tau_vector = 0:0.01:1.0;
meanwait_vector = zeros(length(tau_vector),1);
MinTau_Budget = GetMinTau(S1,S2);
for i = 1 : length(tau_vector)
    meanwait_vector(i) = MeanWait(tau_vector(i),S1,S2);
end
% set(gcf,'color','none'); 
% ylim([0 1])
% xlim([0 1])
plot(tau_vector,meanwait_vector,'r'); 
hold all
plot([0 1],[Epsilon Epsilon],'k','linewidth',2);
plot([MinTau_Budget MinTau_Budget],[0 10],'k','linewidth',2)
% ylim([0 ])
% yt = get(gca, 'YTick');
% set (gca, 'YTickLabel', log(yt)/log(2));
% set(gca, 'YScale', 'log')
% ylim([log(10^(10^(-3))) 10])
plot([MinTau_Budget MinTau_Budget],[0 10],'k','linewidth',2)
ylim([0 max(meanwait_vector)])
% ylim([0 7])
% set(gca,'xtick')
% set(gca,'ytick')
% set(gca,'color','none');
% set(gcf,'InvertHardCopy','off');