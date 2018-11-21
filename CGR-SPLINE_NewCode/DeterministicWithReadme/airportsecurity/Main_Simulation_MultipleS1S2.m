tau = 0:0.005:1;

iseed0 = [  11111111,22222222,33333333,44444444,55555555,...
            66666666,77777777,88888888,99999999,...
            111111111,222222222,3333333333,444444444,555555555,...
            666666666,777777777,888888888,999999999,...
            1111111111,2222222222,3333333333];
%iseed0(1:3) = [1256804093 582283093 2268262506];
num_warmup = 200;
num_simulation = 2000;

for s1 = 1:11
    for s2 = 1:11
        if exist(['E:\AirportModel_Simulation_050616\Results_Model1\New_AvgW_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'file') ~= 2
            continue
        end
        avg_time = zeros(size(tau));
        se_time = zeros(size(tau));
        iseed_vector = iseed0(1:(2+s1+s2));
        for i = 1 : length(tau)
            disp(tau(i))
            [avg_time(i),se_time(i)] = Simulation_AirportModel(tau(i),s1,s2,iseed_vector,num_warmup,num_simulation);
            disp(avg_time(i)*60)
            disp(se_time(i)*60)
        end

        close all
        openfig(['E:\AirportModel_Simulation_050616\Results_Model1\New_AvgW_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'])
        Tmp = get(gca,'children');
        XData = get(Tmp(3),'XData');
        YData = get(Tmp(3),'YData');

        figure
        plot(XData,YData,'k','linewidth',5)
        hold all
        plot(tau,avg_time*60,'b','linewidth',3)
        plot(tau,(avg_time+ 2.576*se_time)*60,'m','linewidth',1)
        plot(tau,(avg_time- 2.576*se_time)*60,'r','linewidth',1)
        title(['Average Waiting Time (S_1 = ' num2str(s1) ', S_2 = ' num2str(s2) ')'])
        ylabel('Average Waiting Time(min)')
        xlabel('tau')
        legend('Mean waiting time','Average Waiting Time (N=50)',' +2.576 se',' -2.576 se')

        saveas(gcf,['Case1_S1_' num2str(s1) '_S2_' num2str(s2) '.png'],'png')
        saveas(gcf,['Case1_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'fig')
        
        figure
        plot(XData,YData,'k','linewidth',5)
        hold all
        plot(tau,avg_time*60,'b','linewidth',3)
        plot(tau,(YData)+ 2.576*se_time*60,'m','linewidth',1)
        plot(tau,(YData)- 2.576*se_time*60,'r','linewidth',1)
        title(['Average Waiting Time (S_1 = ' num2str(s1) ', S_2 = ' num2str(s2) ')'])
        ylabel('Average Waiting Time(min)')
        xlabel('tau')
        legend('Mean waiting time','Average Waiting Time (N=50)',' +2.576 se',' -2.576 se')

        saveas(gcf,['Test_SimulationResults_2\Case1_S1_' num2str(s1) '_S2_' num2str(s2) '.png'],'png')
        saveas(gcf,['Test_SimulationResults_2\Case1_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'fig')
        
    end
end
