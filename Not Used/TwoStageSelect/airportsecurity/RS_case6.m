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
Beta1 = 4;
Beta2 = 38;
Budget = 55;
Epsilon = 8;


global N
global lambda
global mu1
global mu2
global theta
N = 50;
lambda = 20;
mu1 = 9;
mu2 = 6;
theta = 0.25;

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

tautotal = 0;
test = 0;
average = 0;
ba = 5;
totalmobs = 0;
totaltwo = 0;
samplese =0;
totalm= 0;
trial=10;

big_s1  = floor((Budget-Beta1)/c1);
big_s2 = floor((Budget-Beta1)/c2);

tauma = ones(big_s1,big_s2);
sema = ones(big_s1,big_s2);
moma = ones(big_s1,big_s2);
tautoma=zeros(big_s1,big_s2);
testma=zeros(big_s1,big_s2);
tausqrma=zeros(big_s1,big_s2);
trial_best_s = zeros(trial,1);
trial_best_tau = zeros(trial,1);

for t=1:10

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
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(j,k,[111111+i*4+t*100000 222222+i*4+t*100000 333333+i*4+t*100000 444444+i*4+t*100000],...
    2000,2000,...
    1e-3,1e-7)

if tau > 0;
tautotal = tautotal+tau;
totaltwo = totaltwo+tau^2;

test = test+1;
end

totalmobs = totalmobs + mobs;
end

if tautotal >0;

average = tautotal/test;
samplese = ((totaltwo - test*average^2)/(test-1))^0.5;

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

[ snima ] = OCBAf( sencond_tauavma,sencond_sema )


while totalm < 40000000

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
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(j,k,[1111+nose*1000+i*4+t*100000 2222+nose*1000+i*4+t*100000 3333+nose*1000+i*4+t*100000 4444+nose*1000+i*4],...
    200,2000,...
    1e-3,1e-7)

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
 sencond_sema(j+(k-1)*big_s1)=((second_tuasqrma(j+(k-1)*big_s1) - second_test(j+(k-1)*big_s1)*sencond_tauavma(j+(k-1)*big_s1)^2)/(second_test(j+(k-1)*big_s1)-1))^0.5;
 


   end    

sencond_moma(j+(k-1)*big_s1) = totalmobs; 

totalm = totalm + totalmobs;


        end
    end    
end


[ snima ] = OCBAf( sencond_tauavma,sencond_sema )

end

cur_best_tau = 100;
cur_best_s = 0;

for j=1:big_s1
    for k=1:big_s2
        
       if  sencond_tauavma(j+(k-1)*big_s1) < cur_best_tau
           cur_best_tau = sencond_tauavma(j+(k-1)*big_s1);
           cur_best_s = j+(k-1)*big_s1;
           
        end
    end
end

trial_best_s(t,1) = cur_best_s;
trial_best_tau(t,1)= cur_best_tau;


disp(['R&S Case5 result:  BestS = (' 'trial:' num2str(t) ' end, best_s=' num2str(cur_best_s) ])

xlswrite('R&S_Case5_best_s.xls',trial_best_s);
xlswrite('R&S_Case5_best_tau.xls',trial_best_tau);
xlswrite('R&S_Case5_average_tau.xls',sencond_tauavma);
xlswrite('R&S_Case5_SE.xls',sencond_sema);
xlswrite('R&S_Case5_test.xls',second_test);
end


