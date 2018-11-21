close all
cd('C:\Users\BigGhost\Documents\MATLAB\Project_Smoothing\AirportModel')
s1 = 4;
s2 = 4;
load(['C:\Users\BigGhost\Documents\MATLAB\Project_Smoothing\AirportModel\s1_' num2str(s1) '_s2_' num2str(s2) '.mat'])

figure('position',[0,0,600,400])
plot(0.01:0.01:1,WaitingTime_ThisCondition)
title(['Average Waiting Time (S_{1} = ' num2str(s1) ',S_{2} = ' num2str(s2) ')'])
xlabel('tau')
ylabel('Average Waiting Time (Hr)')

saveas(gcf,['AvgW_S1_' num2str(s1) '_S2_' num2str(s2) '.png'],'png')
saveas(gcf,['AvgW_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'fig')


figure('position',[0,0,600,400])
hold all
plot(0.01:0.01:1,Num_NonSelectee_ThisCondition)
plot(0.01:0.01:1,Num_Selectee_ThisCondition)
title(['#Visitors (S_{1} = ' num2str(s1) ',S_{2} = ' num2str(s2) ')'])
xlabel('tau')
ylabel('#Visitors')
legend({'Non-Selectee','Selectee'},'Location','NorthEast')

saveas(gcf,['Visitors_S1_' num2str(s1) '_S2_' num2str(s2) '.png'],'png')
saveas(gcf,['Visirors_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'fig')

figure('position',[0,0,600,400])
hold all
plot(0.01:0.01:1,WaitingTime_NonSelectee_ThisCondition )
plot(0.01:0.01:1,WaitingTime_Selectee_ThisCondition )
title(['Waiting Time (S_{1} = ' num2str(s1) ',S_{2} = ' num2str(s2) ')'])
xlabel('tau')
ylabel('Waiting Time(Hr)')
legend({'Non-Selectee','Selectee'},'Location','NorthEast')

saveas(gcf,['W_S1_' num2str(s1) '_S2_' num2str(s2) '.png'],'png')
saveas(gcf,['W_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'fig')

