
s1 = 2;
s2 = 1;
iseed = 1111111111 ;
num_warmup = 200;
num_simulation = 2000;

testseed = 1 : 100;
avg_time = zeros(size(testseed));
se_time = zeros(size(testseed));

target_tau = 0.12;
for i = 1 : length(testseed)
    disp(i)
    iseed = testseed(i)* ones(1,1+s1+s2);
    [avg_time(i),se_time(i)] = Simulation5_AirportModel(target_tau,s1,s2,iseed,num_warmup,num_simulation);
end

close all
openfig('E:\AirportModel_Simulation_050616\Results_Model1\New_AvgW_S1_2_S2_1.fig')
Tmp = get(gca,'children');

XData = get(Tmp(3),'XData');
YData = get(Tmp(3),'YData');

Idx = find(XData == target_tau);

figure
plot([1 length(testseed)],[YData(Idx) YData(Idx)])
hold all
plot(1:length(testseed),avg_time*60)
legend('Mean Waiting Time','Average Waiting Time')
xlabel('Trial Index')
ylabel('Waiting Time')
title('Waiting Time Comparison (S_1 = 2, S_1 = 1, tau = 0.1)')