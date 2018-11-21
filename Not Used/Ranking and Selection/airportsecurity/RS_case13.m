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
Beta2 = 28;
Budget = 95;
Epsilon = 27;


global N
global lambda
global mu1
global mu2
global theta
N = 50;
lambda = 24;
mu1 = 3;
mu2 = 2;
theta = 0.55;

global tautotal
global test
global average
global ba
global totalmobs
global totaltwo
global samplese
global big_s1
global big_s2
global totalm
global trial
global TotalBudget
global deltaa

TotalBudget = 50000000;
trial=30;
ba = 3;
deltaa = 20;
Tol = 1e-5;
mk = 2000;
warmup = 2000;

tautotal = 0;
test = 0;
average = 0;
totalmobs = 0;
totaltwo = 0;
samplese =0;
totalm= 0;


big_s1  = floor((Budget-Beta1)/c1);
big_s2 = floor((Budget-Beta1)/c2);

tauma = ones(big_s1,big_s2);
sema = ones(big_s1,big_s2);
moma = ones(big_s1,big_s2);
tautoma=zeros(big_s1,big_s2);
testma=zeros(big_s1,big_s2);
tausqrma=zeros(big_s1,big_s2);
trial_best_s = zeros(trial,2);
trial_best_tau = zeros(trial,1);
trial_best_SE = zeros(trial,1); 
trial_best_test = zeros(trial,1);
trial_cost = zeros(trial,1);

for t=1:trial

    totalm=0;
    
for j=1:big_s1
  up_s2 =  floor((Budget-Beta1-c1*j)/c2) ;
  low_s2 = floor(max((Budget-Beta2-c1*j)/c2,1)) ;
  
  
  
    for k=1:big_s2
     
        if k >= low_s2 && k<= up_s2
        
tautotal = 0;
totaltwo = 0;
totalmobs = 0;
average = 0;
test = 0;
samplese=0;


for i=1:ba
      
    
    
    [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(j,k,[111112+i*4+t*100000 222223+i*4+t*100000 333334+i*4+t*100000 444444+i*4+t*100000],...
    warmup,mk,...
    1e-3,Tol)

if tau > 0;
tautotal = tautotal+tau;
totaltwo = totaltwo+tau^2;

test = test+1;
end

totalmobs = totalmobs + mobs;
end

if tautotal >0;

average = tautotal/test;
if test >1
samplese = ((totaltwo - test*average^2)/(test-1))^0.5;
end
    tauma(j+(k-1)*big_s1) = average;
    sema(j+(k-1)*big_s1) = samplese;
    tautoma(j+(k-1)*big_s1) = tautotal;
    tausqrma(j+(k-1)*big_s1) = totaltwo;
    testma(j+(k-1)*big_s1) = test;
end    

moma(j+(k-1)*big_s1) = totalmobs;

totalm = totalm + totalmobs;
        end
    end
end


sencond_tauavma=tauma;
sencond_sema=sema;
sencond_moma=moma;
second_tuasqrma=tausqrma;
second_test=testma;
second_tautoma=tautoma;
nose=0;

[ snima ] = OCBAf( sencond_tauavma,sencond_sema,deltaa)


while totalm < TotalBudget

 nose =  nose+1;
    

for j=1:big_s1
    for k=1:big_s2
        
        
        if snima(j+(k-1)*big_s1)>0
            
            tautotal=0;
            totaltwo=0;
            test=0;
            totalmobs=0;
            
   for i=1:snima(j+(k-1)*big_s1)
      
    
    
    [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(j,k,[1112+nose*1000+i*4+t*100000 2223+nose*1000+i*4+t*100000 3334+nose*1000+i*4+t*100000 4444+nose*1000+i*4],...
    warmup,mk,...
    1e-3,Tol)

          if tau > 0;
               tautotal = tautotal+tau;
               totaltwo = totaltwo+tau^2;
               test = test+1;
          end

         totalmobs = totalmobs + mobs;
   end


  if tautotal >0;

 second_tautoma(j+(k-1)*big_s1) = tautotal+second_tautoma(j+(k-1)*big_s1);
 second_test(j+(k-1)*big_s1)=test+second_test(j+(k-1)*big_s1);
 sencond_tauavma(j+(k-1)*big_s1) = second_tautoma(j+(k-1)*big_s1)/second_test(j+(k-1)*big_s1);
 
 
 second_tuasqrma(j+(k-1)*big_s1)=second_tuasqrma(j+(k-1)*big_s1)+totaltwo;
 
 if second_test(j+(k-1)*big_s1) >1
 sencond_sema(j+(k-1)*big_s1)=((second_tuasqrma(j+(k-1)*big_s1) - second_test(j+(k-1)*big_s1)*sencond_tauavma(j+(k-1)*big_s1)^2)/(second_test(j+(k-1)*big_s1)-1))^0.5;
 end


   end    

sencond_moma(j+(k-1)*big_s1) = totalmobs; 

totalm = totalm + totalmobs;


        end
    end    
end


[ snima ] = OCBAf( sencond_tauavma,sencond_sema,deltaa)

end

cur_best_tau = 100;
cur_best_s = 0;
cur_best_j = 0;
cur_best_k = 0;

for j=1:big_s1
    for k=1:big_s2
        
       if  sencond_tauavma(j+(k-1)*big_s1) < cur_best_tau
           cur_best_tau = sencond_tauavma(j+(k-1)*big_s1);
           cur_best_s = j+(k-1)*big_s1;
           cur_best_j = j;
           cur_best_k = k;
        end
    end
end

 [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(cur_best_j,cur_best_k,[11111111 22222222 33333333 44444444],...
    200,48600,...
    1e-3,0.005/(48600^(1/2)))

trial_best_s(t,1) = cur_best_j;
trial_best_s(t,2) = cur_best_k;
if tau > 0
trial_best_tau(t,1)= tau;
end
trial_best_SE(t,1) = sencond_sema(cur_best_s);
trial_best_test(t,1) = second_test(cur_best_s);
trial_cost(t,1) = totalm+ mobs;

disp(['R&S Case13 result:  BestS = (' 'trial:' num2str(t) ' end, best_s=' num2str(cur_best_s) ])

xlswrite('R&S_Case13_best_s.xls',trial_best_s);
xlswrite('R&S_Case13_best_tau.xls',trial_best_tau);
xlswrite('R&S_Case13_best_SE.xls',trial_best_SE);
xlswrite('R&S_Case13_best_test.xls',trial_best_test);
xlswrite('R&S_Case13_TotalCost.xls',trial_cost);
xlswrite(strcat('R&S_Case13_average_tau_Trial',int2str(t),'.xls'),sencond_tauavma);
xlswrite(strcat('R&S_Case13_SE_Trial',int2str(t),'.xls'),sencond_sema);
xlswrite(strcat('R&S_Case13_test_Trial',int2str(t),'.xls'),second_test);
end


