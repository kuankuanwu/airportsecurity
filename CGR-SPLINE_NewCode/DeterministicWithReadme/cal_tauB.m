% global c1
% global c2
% global Beta1
% global Beta2
% global Budget
% global theta
% global N
% global Epsilon
% Epsilon = 25;
% N=200;
% c1=2;
% c2=5;
% Beta1=5;
% Beta2=50;
% Budget=95;
% theta = 0.55;
% global mu1
% global mu2
% global Lambda
% global RefMatrix
% RefMatrix=-ones(100,100);  
% mu1=0.66667;
% mu2=0.016667;
% Lambda=2;
% Range_S1 = floor((Budget-Beta1)/c1);
% Range_S2 = floor((Budget-Beta1)/c2);
% 
% TauB = -ones(Range_S1,Range_S2);
% time_matrix=-ones(Range_S1,Range_S2);
% 
% for s1 = 0 : Range_S1
%     UP_S1 = floor((Budget-Beta1-c1*s1)/c2);
%     LO_S1=  max(floor((Budget-Beta2-c1*s1)/c2),0);
%     
%     for s2 = LO_S1 :UP_S1
%         if  s2 == 0
%             continue
%         end
% %         mintau = GetMinTau(s1,s2);
% %         
% %         TauB(s1,s2) = mintau;
%             [tau,MeanWaitingTime,MinTau_Budget,flag] = GetOptimTau(s1,s2,1);
%             if isempty(tau)==1
%                 continue
%             end
%             TauB(s1,s2) = tau;
% %         MeanWaitingTime=Deterministic_TestBlockTridiagonal_LargeMatrix2(mintau,s1,s2);
% %         time_matrix(s1,s2)=MeanWaitingTime;
%         
%     end
%     
% end

global c1
global c2
global Beta1
global Beta2
global Budget
global theta
c1=2;
c2=3;
Beta1=5;
Beta2=35;
Budget=95;
theta = 0.55;
Range_S1 = floor((Budget-Beta1)/c1);
Range_S2 = floor((Budget-Beta1)/c2);

TauB = -ones(Range_S1,Range_S2);
for s1 = 1 : Range_S1
    UP_S1 = floor((Budget-Beta1-c1*s1)/c2);
    LO_S1=  max(floor((Budget-Beta2-c1*s1)/c2),0);
    
    for s2 = LO_S1 :UP_S1
        if  s2 == 0
            continue
        end
        mintau = GetMinTau(s1,s2);
        
        TauB(s1,s2) = mintau;
    end
    
end