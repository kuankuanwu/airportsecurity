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
Epsilon = 8.8;


global N
global lambda
global mu1
global mu2
global theta
N = 50;
lambda = 24;
mu1 = 8;
mu2 = 4;
theta = 0.0625;

S1 = 9;
S2 = 5;

mk=1000;
warm=100;

optimaltau = zeros(20,50);
mobsA = zeros(20,50);

for i=1:1000
   
    [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(S1,S2,[1112+i*4 2223+i*4 3334+i*4 4444+i*4],...
    warm,mk,...
    1e-3,1e-5);
   
optimaltau(ceil(i/50),(i-(ceil(i/50)-1)*50)) = tau;
mobsA(ceil(i/50),(i-(ceil(i/50)-1)*50)) =mobs;

end
save(strcat('Case5.2_mobs_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'mobsA')
save(strcat('Case5.2_optimal_tau_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'optimaltau')


mk=1000;
warm=1000;

optimaltau = zeros(20,50);
mobsA = zeros(20,50);

for i=1:1000
   
    [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(S1,S2,[1112+i*4 2223+i*4 3334+i*4 4444+i*4],...
    warm,mk,...
    1e-3,1e-5);
   
optimaltau(ceil(i/50),(i-(ceil(i/50)-1)*50)) = tau;
mobsA(ceil(i/50),(i-(ceil(i/50)-1)*50)) =mobs;
end

save(strcat('Case5.2_mobs_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'mobsA')
save(strcat('Case5.2_optimal_tau_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'optimaltau')

mk=2000;
warm=100;

optimaltau = zeros(20,50);
mobsA = zeros(20,50);

for i=1:1000
   
    [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(S1,S2,[1112+i*4 2223+i*4 3334+i*4 4444+i*4],...
    warm,mk,...
    1e-3,1e-5);
   
optimaltau(ceil(i/50),(i-(ceil(i/50)-1)*50)) = tau;
mobsA(ceil(i/50),(i-(ceil(i/50)-1)*50)) =mobs;

end

save(strcat('Case5.2_mobs_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'mobsA')
save(strcat('Case5.2_optimal_tau_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'optimaltau')

mk=2000;
warm=1000;

optimaltau = zeros(20,50);
mobsA = zeros(20,50);

for i=1:1000
   
    [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(S1,S2,[1112+i*4 2223+i*4 3334+i*4 4444+i*4],...
    warm,mk,...
    1e-3,1e-5);
   
optimaltau(ceil(i/50),(i-(ceil(i/50)-1)*50)) = tau;
mobsA(ceil(i/50),(i-(ceil(i/50)-1)*50)) =mobs;

end

save(strcat('Case5.2_mobs_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'mobsA')
save(strcat('Case5.2_optimal_tau_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'optimaltau')

mk=5000;
warm=100;

optimaltau = zeros(20,50);
mobsA = zeros(20,50);

for i=1:1000
   
    [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(S1,S2,[1112+i*4 2223+i*4 3334+i*4 4444+i*4],...
    warm,mk,...
    1e-3,1e-5);
   
optimaltau(ceil(i/50),(i-(ceil(i/50)-1)*50)) = tau;
mobsA(ceil(i/50),(i-(ceil(i/50)-1)*50)) =mobs;

end

save(strcat('Case5.2_mobs_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'mobsA')
save(strcat('Case5.2_optimal_tau_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'optimaltau')

mk=5000;
warm=1000;

optimaltau = zeros(20,50);
mobsA = zeros(20,50);

for i=1:1000
   
    [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(S1,S2,[1112+i*4 2223+i*4 3334+i*4 4444+i*4],...
    warm,mk,...
    1e-3,1e-5);
   
optimaltau(ceil(i/50),(i-(ceil(i/50)-1)*50)) = tau;
mobsA(ceil(i/50),(i-(ceil(i/50)-1)*50)) =mobs;

end

save(strcat('Case5.2_mobs_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'mobsA')
save(strcat('Case5.2_optimal_tau_mk = ',int2str(mk),'warmup = ',int2str(warm),'.mat'),'optimaltau')

