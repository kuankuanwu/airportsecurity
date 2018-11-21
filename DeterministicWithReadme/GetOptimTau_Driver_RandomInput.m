
InputFileName = 'input.txt';
outfile = fopen('Output.txt','wt');
InitializeAndReadInput(InputFileName)


%(Constraint)
global c1
global c2
global Beta1
global Beta2
global Budget
global Epsilon


%(Waiting Time)
global N
global Lambda
global mu1
global mu2
global theta


Range_S1 = floor((Budget-Beta1)/c1);
Range_S2 = floor((Budget-Beta1)/c2);

OptimTau = -ones(Range_S1,Range_S2);
OptimTau_S1Zero = -ones(1,Range_S2);
OptimTau_S2Zero = -ones(Range_S1,1);

OptimTau_MinMeanWaitingTime = -ones(Range_S1,Range_S2);
OptimTau_MinMeanWaitingTime_S1Zero = -ones(1,Range_S2);
OptimTau_MinMeanWaitingTime_S2Zero = -ones(Range_S1,1);

trial=0;
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

while 1 == 1
while 1 ==1
 Lambda = Lambda + randi([-10 10]);
 mu1 = mu1 + randi([-5 5]);
 mu2 = mu2 + randi([-5 5]);   
 Epsilon = Epsilon +  randi([-10 10]); 

 if Lambda*mu1*mu2* Epsilon > 0
    break
 end   
 
end    
trial=trial+1;
filename = ['Trial' num2str(trial) '.mat'];

for s1 = 0: Range_S1
    UP_S1 = floor((Budget-Beta1-c1*s1)/c2);
    LO_S1=  max(floor((Budget-Beta2-c1*s1)/c2),0);
    for s2 = LO_S1: UP_S1
        if s1 == 0 && s2 == 0
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
        save(filename,'OptimTau_S1Zero','OptimTau_MinMeanWaitingTime_S1Zero',...
            'OptimTau_S2Zero','OptimTau_MinMeanWaitingTime_S2Zero',...
            'OptimTau','OptimTau_MinMeanWaitingTime');
        %xlswrite('ALL_optimaltau_Case9.xls',OptimTau);
        
    end
end

fprintf(outfile,'Trial (%d)\n',trial );
fprintf(outfile,'Lambda = (%d)\n', Lambda );
fprintf(outfile,'mu1 = (%d)\n', mu1 );
fprintf(outfile,'mu2 = (%d)\n', mu2 );
fprintf(outfile,'Epsilon = (%d)\n', Epsilon );

besttau=1;
bests1 = 0;
bests2=0;
for s1 = 1: Range_S1
for s2 = LO_S1: UP_S1
    if  OptimTau(s1,s2) <besttau
    if s2 >0
   besttau =  OptimTau(s1,s2);
   bests1=s1;
   bests2=s2;
    end    
    end
end
end
fprintf(outfile,'BestTau = (%d)\n', besttau );
for s1 = 1: Range_S1-1
for s2 = LO_S1: UP_S1
     
    if OptimTau(s1,s2)~=-1
    if s2>0   
        if OptimTau(s1,s2) < OptimTau(s1+1,s2)
        if OptimTau(s1,s2) < OptimTau(s1-1,s2)
        if OptimTau(s1,s2) < OptimTau(s1,s2+1)
        if OptimTau(s1,s2) < OptimTau(s1,s2-1)
            fprintf(outfile,'Local optimal - Global optimal  (%d)\n',OptimTau(s1,s2)-besttau );   
        end
        end
        end
        end
    end    
    end 
end
end

end
