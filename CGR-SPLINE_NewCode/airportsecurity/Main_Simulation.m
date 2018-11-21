tau = 0:0.005:1;
s1 = 1;
s2 = 1;
iseed = 2222222222;
num_warmup = 200;
num_simulation = 2000;


avg_time = zeros(size(tau));
se_time = zeros(size(tau));
Length1 = zeros(size(tau));
Length2 = zeros(size(tau));
for i = 1 : length(tau)
    disp(tau(i))
    iseed_vector = iseed*ones(1,1+s1+s2);
    [avg_time(i),se_time(i)] = Simulation6b_AirportModel(tau(i),s1,s2,iseed_vector,num_warmup,num_simulation);
end

close all
openfig('E:\AirportModel_Simulation_050616\Results_Model1\New_AvgW_S1_2_S2_1.fig')
Tmp = get(gca,'children');
XData = get(Tmp(3),'XData');
YData = get(Tmp(3),'YData');

figure
plot(tau,avg_time*60,'b','linewidth',5)
hold all
plot(XData,YData,'k','linewidth',3)
plot(tau,(avg_time-2.576*se_time)*60,'r','linewidth',1)
plot(tau,(avg_time+2.576*se_time)*60,'g','linewidth',1)
title(['Average Waiting Time (S_1 = ' num2str(s1) ', S_2 = ' num2str(s2) ', N = 50)'])
ylabel('Average Waiting Time(min)')
xlabel('tau')

% figure
% plot(tau,Length1)
% hold all
% plot(tau,Length2)
% title(['Average Queue Length (S_1 = ' num2str(s1) ', S_2 = ' num2str(s2) ', N = 50)'])
% ylabel('Average Queue Length')
% xlabel('tau')
