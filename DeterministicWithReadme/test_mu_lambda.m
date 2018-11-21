global c1
global c2
global Beta1
global Beta2
global Budget
global theta
c1=1;
c2=1;
Beta1=5;
Beta2=25;
Budget=100;
theta = 0.55;
global mu1
global mu2
global Lambda 
global Epsilon
Epsilon=12/60;
Lambda=40;
mu1=8;
mu2=4;
global N
N=50;
Range_S1 = floor((Budget-Beta1)/c1);
Range_S2 = floor((Budget-Beta1)/c2);
global RefMatrix
    RefMatrix = -ones(Range_S1,Range_S2); 

test_tau =-ones(Range_S1,Range_S2);
test_wait=-ones(Range_S1,Range_S2);
for s1 = 18:37   
    UP_S1 = floor((Budget-Beta1-c1*s1)/c2);
    LO_S1=  max(floor((Budget-Beta2-c1*s1)/c2),0);
    for s2 = LO_S1 :UP_S1
        [tau,MeanWaitingTime,MinTau_Budget,flag] = GetOptimTau(s1,s2,1);
        if flag > 0
        test_tau(s1,s2)=tau;
        test_wait(s1,s2)=MeanWaitingTime;
    
        end
        
    end
end
    