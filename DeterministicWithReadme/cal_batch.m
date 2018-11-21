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
 S1=4;
 S2=22;
 tau=0.281730242;
 iseed =  [1653381726,2882905302,1955172277,44444444];
 Num_Simulation = 54000 ;
% 
h1=[];
h2=[];
res=[];

batch_vector = 20:10:1000;
se_vector = zero s(length(batch_vector),1);
for i = 1 : length(tau_vector)
    [MeanWaitingTimeSim,SE]=Simulation_AirportModel(tau,S1,S2,iseed,Num_Warmup,Num_Simulation);
    se_vector(i) = SE;
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
set(gca, 'YScale', 'log')
% ylim([log(10^(10^(-5)))  10 ])
plot([MinTau_Budget MinTau_Budget],[0 10],'k','linewidth',2)
% set(gcf, 'Position', [100, 100, 100, 100])
pbaspect([1 1 1])
% ylim([0 max(meanwait_vector)])
% ylim([0 7])
% set(gca,'xtick')
% set(gca,'ytick')
% set(gca,'color','none');
% set(gcf,'InvertHardCopy','off');