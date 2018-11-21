%(Constraint)
global c1
global c2
global Beta1
global Beta2
global Budget
global Epsilon
c1 = 3;
c2 = 5;
Beta1 = 5;
Beta2 = 20;
Budget = 44;
Epsilon = 1;

%(Waiting Time)
global N
global Lambda
global mu1
global mu2
global theta
N = 50;
Lambda = 18;
mu1 = 8;
mu2 = 6;
theta = 0.0625;


Range = 1:10;

FinalMinMeanWaitingTime = -ones(11,11);
FinalMinMeanWaitingTime_S1Zero = -ones(11,11);
FinalMinMeanWaitingTime_S2Zero = -ones(11,11);

FinalTau_MinMeanWaitingTime = -ones(11,11);
FinalTau_MinMeanWaitingTime_S1Zero = -ones(11,11);
FinalTau_MinMeanWaitingTime_S2Zero = -ones(11,11);

for s1 = 1: 1
    for s2 = 1 :1
        if s1 == 0 && s2 == 0
            continue
        end
        close all
        if s1 == 0
            tau = 0;
        elseif s2 == 0
            tau = 1;
        else
            tau = 0:0.005:1;
        end
        MinTau_Budget = GetMinTau(s1,s2);
        if ~isreal(MinTau_Budget)
            continue
        end
        if MinTau_Budget < 0
            MinTau_Budget = 0;
        end
        WaitingTime_ThisCondition = zeros(size(tau));
        PE1 = zeros(size(tau));
        PE2 = zeros(size(tau));
        PE12 = zeros(size(tau));
        p = zeros(size(tau));
        for i = 1 : length(tau)
            disp(tau(i))
            p(i) = (-exp(-tau(i)/theta) + 1)/(1 - exp(-1/theta));
            [WaitingTime_ThisCondition(i),PE1(i),PE2(i),PE12(i)] = MeanWait(tau(i),s1,s2);   
        end
        
        if s1 == 0
            FinalTau_MinMeanWaitingTime_S1Zero(s2) = 0;
            FinalMinMeanWaitingTime_S1Zero(s2) = WaitingTime_ThisCondition;
        elseif s2 == 0
            FinalTau_MinMeanWaitingTime_S2Zero(s1) = 1;
            FinalMinMeanWaitingTime_S2Zero(s1) = WaitingTime_ThisCondition;
        else
            Idx = find(WaitingTime_ThisCondition == min(WaitingTime_ThisCondition));
            Idx = Idx(1);
            FinalTau_MinMeanWaitingTime(s1,s2) = tau(Idx);
            FinalMinMeanWaitingTime(s1,s2) = WaitingTime_ThisCondition(Idx);
        end
        %MinTau_Budget = 0.5;
        figure('position',[0,0,600,400])
        if length(tau) == 1
            plot(tau,WaitingTime_ThisCondition,'ro-')
        else
            plot(tau,WaitingTime_ThisCondition,'r-')
        end
        
        hold all
        plot([0 1],[Epsilon Epsilon],'k','linewidth',3);
        plot([MinTau_Budget MinTau_Budget],[0 max(Epsilon*2,max(WaitingTime_ThisCondition))],'k','linewidth',5)
        xlim([0 1]);
        title(['Mean Waiting Time (S_{1} = ' num2str(s1) ',S_{2} = ' num2str(s2) ')'])
        xlabel('tau')
        ylabel('Mean Waiting Time (Min)')

        saveas(gcf,['AvgW_S1_' num2str(s1) '_S2_' num2str(s2) '.png'],'png')
        saveas(gcf,['AvgW_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'fig')
        
        figure('position',[0,0,600,400])
        if length(tau) == 1
            plot(tau,p.*(PE1 + PE12),'ro-')
            hold all
            plot(tau,(1-p).*(PE2 + PE12),'bo-')
        else
            plot(tau,p.*(PE1 + PE12),'r-')
            hold all
            plot(tau,(1-p).*(PE2 + PE12),'b-')
        end
        xlim([0 1]);
        title(['Enterable (S_{1} = ' num2str(s1) ',S_{2} = ' num2str(s2) ')'])
        legend('P(Non-Selected:Enterable)','P(Selected:Enterable)')
        xlabel('tau')
        ylabel('P')

        saveas(gcf,['Enterable_S1_' num2str(s1) '_S2_' num2str(s2) '.png'],'png')
        saveas(gcf,['Enterable_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'fig')
        
        p = (-exp(-tau/theta) + 1)/(1 - exp(-1/theta));
        TotalLambda = Lambda*p.*(PE1 + PE12) + Lambda*(1-p).*(PE2 + PE12);
        
        figure('position',[0,0,600,400])
        if length(tau) == 1
            plot(tau,TotalLambda,'ro-')
        else
            plot(tau,TotalLambda,'r-')
        end
        xlim([0 1]);
        title(['TotalArrivalRate: (S_{1} = ' num2str(s1) ',S_{2} = ' num2str(s2) ')'])
        legend('Arrival Rate')
        xlabel('tau')
        ylabel('Arrival Rate')

        saveas(gcf,['ArrivalRate_S1_' num2str(s1) '_S2_' num2str(s2) '.png'],'png')
        saveas(gcf,['ArrivalRate_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'fig')
        
    end
end


