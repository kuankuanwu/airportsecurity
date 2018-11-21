tau = 0:0.005:1;

iseed0 = [  11111111,22222222,33333333,44444444,55555555,...
            66666666,77777777,88888888,99999999,...
            111111111,222222222,3333333333,444444444,555555555,...
            666666666,777777777,888888888,999999999,...
            1111111111,2222222222,3333333333];
num_warmup = 2000;
num_simulation = 10000;

d1 = 0.7;
d2 = 0.98;

for s1 = 1:1
    for s2 = 1:1
        if exist(['E:\AirportModel_Simulation_050616\Results_Model1\New_AvgW_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'file') ~= 2
            continue
        end
        avg_time = zeros(size(tau));
        se_time = zeros(size(tau));
        iseed_vector = iseed0(1:(2+s1+s2));
        AvgR1 = zeros(size(tau));
        AvgR2 = zeros(size(tau));
        Avg_SL = zeros(size(tau));
        Se_SL = zeros(size(tau));
        EvalR1 = zeros(size(tau));
        EvalR2 = zeros(size(tau));
        Eval_SL = zeros(size(tau));
        for i = 1 : length(tau)
            R1_Sim = zeros(1,20);
            R2_Sim = zeros(1,20);
            SL_Sim = zeros(1,20);
            for j = 1 : 20
                [avg_time(i),se_time(i),Dummy,R1_Sim(j),R2_Sim(j)] = Simulation_AirportModel(tau(i),s1,s2,iseed_vector,num_warmup,num_simulation);
                SL_Sim(j) = R1_Sim(j)*d1 + R2_Sim(j)*d2;
            end
            AvgR1(i) = mean(R1_Sim);
            AvgR2(i) = mean(R2_Sim);
            Avg_SL(i) = mean(SL_Sim);
            Se_SL(i) = std(SL_Sim)/sqrt(20);
            

            EvalR1(i) = CalcR1(tau(i),0.0625);
            EvalR2(i) = CalcR2(tau(i),0.0625);
            Eval_SL(i) = d1*EvalR1(i) + d2*EvalR2(i);
            
            disp(tau(i));
            disp(Avg_SL(i));
            disp(Eval_SL(i));
        end

        close all
        
        figure
        plot(tau,Eval_SL,'k','linewidth',5)
        hold all
        plot(tau,Avg_SL,'b','linewidth',3)
        plot(tau,(Avg_SL+ 2.576*Se_SL),'m','linewidth',1)
        plot(tau,(Avg_SL- 2.576*Se_SL),'r','linewidth',1)
        title(['Security Level (S_1 = ' num2str(s1) ', S_2 = ' num2str(s2) ')'])
        ylabel('Security Level')
        xlabel('tau')
        legend('Ori Security Level','Average Security Level',' +2.576 se',' -2.576 se')

        saveas(gcf,['SecurityLevel_S1_' num2str(s1) '_S2_' num2str(s2) '.png'],'png')
        saveas(gcf,['SecurityLevel_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'fig')
        
    end
end
