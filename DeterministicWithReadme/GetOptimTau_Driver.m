for i=1:1
InputFileName = ['input',num2str(i),'.txt'];
InitializeAndReadInput(InputFileName)

%(Constraint)
global c1
global c2
global Beta1
global Beta2
global Budget
global Epsilon
global Tol 
global MaxIter 
Tol = 1E-7;
MaxIter = 1000;

%(Waiting Time)
global N
global Lambda
global mu1
global mu2
global theta
global chainflag
chainflag=1;

Range_S1 = floor((Budget-Beta1)/c1);
Range_S2 = floor((Budget-Beta1)/c2);

OptimTau = -ones(Range_S1,Range_S2);
OptimTau_S1Zero = -ones(1,Range_S2);
OptimTau_S2Zero = -ones(Range_S1,1);

OptimTau_MinMeanWaitingTime = -ones(Range_S1,Range_S2);
OptimTau_MinMeanWaitingTime_S1Zero = -ones(1,Range_S2);
OptimTau_MinMeanWaitingTime_S2Zero = -ones(Range_S1,1);


global RefMatrix
    RefMatrix = -ones(Range_S1,Range_S2); 
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

FileName = ['input',num2str(i),'.txt'];

FileName = ['Case33_N200_New111_' num2str(i) '_c1_' ,num2str(c1),'_c2_',num2str(c2),'_beta1_',num2str(Beta1),'_beta2_',num2str(Beta2),'_Budget_',num2str(Budget),'_epsilon_',num2str(Epsilon),...
    '_lambda_',num2str(Lambda),'_mu1_',num2str(mu1),'_mu2_',num2str(mu2),'.mat']
disp(FileName)

for s1 = 0: Range_S1
    UP_S1 = floor((Budget-Beta1-c1*s1)/c2);
    LO_S1=  max(floor((Budget-Beta2-c1*s1)/c2),0);
    for s2 = LO_S1: UP_S1
        if s1 == 20 && s2 >= 15
            continue
        end
%         if s1 < Prev_Range_S1 && s2 < Prev_Range_S2 && (s1 > 0 && s2 > 0)
%             continue
%         end
        disp(['[' num2str(s1) ',' num2str(s2) ']'])
        [tau, MeanWaitingTime, MinTau_Budget, flag] = GetOptimTau(s1,s2,1);
        if flag > 0
            if s1 == 0
                OptimTau_S1Zero(s2) = tau;
                OptimTau_MinMeanWaitingTime_S1Zero(s2) = MeanWaitingTime;
                disp(['[' num2str(s1) ',' num2str(s2) ']end,optimal tau=' num2str(tau) ])
            elseif s2 == 0
                OptimTau_S2Zero(s1) = tau;
                OptimTau_MinMeanWaitingTime_S2Zero(s1) = MeanWaitingTime;
                disp(['[' num2str(s1) ',' num2str(s2) ']end,optimal tau=' num2str(tau) ])
            else
                OptimTau(s1,s2) = tau;
                OptimTau_MinMeanWaitingTime(s1,s2) = MeanWaitingTime;
                disp(['[' num2str(s1) ',' num2str(s2) ']end,optimal tau=' num2str(tau) ])
            end
        end
        save(FileName,'OptimTau_S1Zero','OptimTau_MinMeanWaitingTime_S1Zero',...
            'OptimTau_S2Zero','OptimTau_MinMeanWaitingTime_S2Zero',...
            'OptimTau','OptimTau_MinMeanWaitingTime');
%         xlswrite('ALL_optimaltau_Case18.xls',OptimTau);
        
    end
end
end
