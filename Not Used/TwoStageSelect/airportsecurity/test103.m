global log_test
log_test_name='Debug_log_SimpleBisection_Test_Case1_RandomInitial.txt';
log_test=fopen(log_test_name, 'w'); 
global c1
global c2
global Beta1
global Beta2
global Budget
global Epsilon
c1 = 2;
c2 = 3;
Beta1 = 5;
Beta2 = 35;
Budget = 44;
Epsilon = 9;


global N
global lambda
global mu1
global mu2
global theta
N = 50;
lambda = 36;
mu1 = 12;
mu2 = 6;
theta = 0.25;

tauma=zeros(30,1);

for i = 1:30

[tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(7,2,[1111+i,2222+i 3333+i 4444+i],...
    200,2000,...
    1e-3,1e-7)

tauma(i,1)=tau;

disp(['i=' num2str(i) ',tau=' num2str(tau)])
xlswrite('TEST.xls',tauma);
end