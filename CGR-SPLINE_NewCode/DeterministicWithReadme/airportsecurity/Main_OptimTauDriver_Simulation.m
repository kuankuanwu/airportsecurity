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
Epsilon = 13;

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

iseed0 = [  11111111,22222222,33333333,44444444,55555555,...
            66666666,77777777,88888888,99999999,...
            111111111,222222222,3333333333,444444444,555555555,...
            666666666,777777777,888888888,999999999,...
            1111111111,2222222222,3333333333];
num_warmup = 200;
num_simulation = 2000;
maxlayer_search = 7;
maxlayer_solve = 5;
tol = 1e-4;

%Range_S1 = 3;
%Range_S2 = 3;
Range_S1 = ceil(Budget/c1);
Range_S2 = ceil(Budget/c2);

OptimTau = -ones(Range_S1,Range_S2);
OptimTau_S1Zero = -ones(1,Range_S2);
OptimTau_S2Zero = -ones(Range_S1,1);

OptimTau_MinMeanWaitingTime = -ones(Range_S1,Range_S2);
OptimTau_MinMeanWaitingTime_S1Zero = -ones(1,Range_S2);
OptimTau_MinMeanWaitingTime_S2Zero = -ones(Range_S1,1);

% Prev_Range_S1 = 11;
% Prev_Range_S2 = 11;
% OptimTau_Tmp = -ones(Range_S1,Range_S2);
% OptimTau_S1Zero_Tmp = -ones(Range_S1,Range_S2);
% OptimTau_S2Zero_Tmp = -ones(Range_S1,Range_S2);
% 
% OptimTau_MinMeanWaitingTime_Tmp = -ones(Range_S1,Range_S2);
% OptimTau_MinMeanWaitingTime_S1Zero_Tmp = -ones(1,Range_S2);
% OptimTau_MinMeanWaitingTime_S2Zero_Tmp = -ones(Range_S1,1);
% 
% OptimTau_Tmp(1:Prev_Range_S1,1:Prev_Range_S2) = OptimTau;
% %OptimTau_S1Zero_Tmp(1:Prev_Range_S2) = OptimTau_S1Zero;
% %OptimTau_S2Zero_Tmp(1:Prev_Range_S1) = OptimTau_S2Zero;
% 
% OptimTau_MinMeanWaitingTime_Tmp(1:Prev_Range_S1,1:Prev_Range_S2) = OptimTau_MinMeanWaitingTime;
% OptimTau_MinMeanWaitingTime_S1Zero_Tmp(1:Prev_Range_S2) = OptimTau_MinMeanWaitingTime_S1Zero;
% OptimTau_MinMeanWaitingTime_S2Zero_Tmp(1:Prev_Range_S1) = OptimTau_MinMeanWaitingTime_S2Zero;
% 
% OptimTau = OptimTau_Tmp;
% OptimTau_S1Zero = OptimTau_S1Zero_Tmp;
% OptimTau_S2Zero = OptimTau_S2Zero_Tmp;
% OptimTau_MinMeanWaitingTime = OptimTau_MinMeanWaitingTime_Tmp;
% OptimTau_MinMeanWaitingTime_S1Zero = OptimTau_MinMeanWaitingTime_S1Zero_Tmp;
% OptimTau_MinMeanWaitingTime_S2Zero = OptimTau_MinMeanWaitingTime_S2Zero_Tmp;




for s1 = 0 : Range_S1
    for s2 = 0 : Range_S2
%for s1 = 3:3
%    for s2 = 3:3
        if s1 == 0 && s2 == 0
            continue
        end
%         if s1 < Prev_Range_S1 && s2 < Prev_Range_S2 && (s1 > 0 && s2 > 0)
%             continue
%         end
        disp(['[' num2str(s1) ',' num2str(s2) ']'])
        [tau, MeanWaitingTime, MinTau_Budget, flag] = GetOptimTau_Simulation(s1,s2,iseed0,num_warmup,num_simulation,maxlayer_search,maxlayer_solve,tol);
        if flag > 0
            if s1 == 0
                OptimTau_S1Zero(s2) = tau;
                OptimTau_MinMeanWaitingTime_S1Zero(s2) = MeanWaitingTime;
            elseif s2 == 0
                OptimTau_S2Zero(s1) = tau;
                OptimTau_MinMeanWaitingTime_S2Zero(s1) = MeanWaitingTime;
            else
                OptimTau(s1,s2) = tau;
                OptimTau_MinMeanWaitingTime(s1,s2) = MeanWaitingTime;
            end
        end
        save('Tmp.mat','OptimTau_S1Zero','OptimTau_MinMeanWaitingTime_S1Zero',...
            'OptimTau_S2Zero','OptimTau_MinMeanWaitingTime_S2Zero',...
            'OptimTau','OptimTau_MinMeanWaitingTime');
    end
end


for s1 = 0 : Range_S1
    for s2 = 0 : Range_S2
%for s1 = 3:3
%    for s2 = 3:3
        if s1 == 0 && s2 == 0
            continue
        end
%         if s1 < Prev_Range_S1 && s2 < Prev_Range_S2 && (s1 > 0 && s2 > 0)
%             continue
%         end
        disp(['[' num2str(s1) ',' num2str(s2) ']'])
        if exist(['E:\AirportModel_Simulation_CGRSpline_062816\SimulationResult_1\New2_DifferentSeed_AvgWaitingTime_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'file') ~= 2
            continue
        end
        openfig(['E:\AirportModel_Simulation_CGRSpline_062816\SimulationResult_1\New2_DifferentSeed_AvgWaitingTime_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'])
        hold all
        plot([0 1],[Epsilon Epsilon],'r-','linewidth',2);
        if s1 == 0
            CurOptimTau = OptimTau_S1Zero(s2);
            CurAvgWT = OptimTau_MinMeanWaitingTime_S1Zero(s2);
        elseif s2 == 0
            CurOptimTau = OptimTau_S2Zero(s1);
            CurAvgWT = OptimTau_MinMeanWaitingTime_S2Zero(s1);
        else
            CurOptimTau = OptimTau(s1,s2);
            CurAvgWT = OptimTau_MinMeanWaitingTime(s1,s2);
        end
        if CurOptimTau > 0
            plot(CurOptimTau,CurAvgWT,'go','markerface','g','markersize',6)
            title(['Averag Waiting Time (S_1 = ' num2str(s1) ', S_2 = ' num2str(s2) '), OptimTau = ' num2str(CurOptimTau)])
        else
            title(['Averag Waiting Time (S_1 = ' num2str(s1) ', S_2 = ' num2str(s2) '), Not Feasible'])
        end
        saveas(gcf,['E:\AirportModel_Simulation_CGRSpline_062816\OptimTau_Simulation_S1_' num2str(s1) '_S2_' num2str(s2) '.png'],'png')
        saveas(gcf,['E:\AirportModel_Simulation_CGRSpline_062816\OptimTau_Simulation_S1_' num2str(s1) '_S2_' num2str(s2) '.fig'],'fig')
    end
end
