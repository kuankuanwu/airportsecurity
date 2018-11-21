global log_test
log_test_name='Debug_log_SimpleBisection_Test_Case1_RandomInitial.txt';
log_test=fopen(log_test_name, 'w'); 
totalbudget_test_name='bugdettest.txt';
buget_test=fopen(totalbudget_test_name,'w');
global c1
global c2
global Beta1
global Beta2
global Budget
global Epsilon
c1 = 1;
c2 = 2;
Beta1 = 5;
Beta2 = 19;
Budget = 62;
Epsilon = 31;


global N
global lambda
global mu1
global mu2
global theta
N = 50;
lambda = 20;
mu1 = 3;
mu2 = 1.9;
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
global TotalBudget
global deltaa


global flag
global changenumber 
global times
times = 2;
flag = 1;
changenumber = 0;
TotalBudget = 50000000;
trial=30;
ba = 3;
deltaa = 20; 
Tol = 1e-5;
m_0 = 200;
warmup = 1000;

global alpha
system_number = 408;
sample = 20;
alpha = 0.05;



tautotal = 0;
test = 0;
average = 0;
totalmobs = 0;
totaltwo = 0;
samplese =0;
totalm= 0;


big_s1  = floor((Budget-Beta1)/c1);
big_s2 = floor((Budget-Beta1)/c2);
testnumber=ones(big_s1,big_s2);
stage_select=[];
stage_output=zeros(50,4);
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
trial_total = zeros(trial,10);
star_trial=4;
end_trial=4;
for t=star_trial:end_trial
  
    totalm=0;
    stage=1;
    for j=1:big_s1
        up_s2 =  floor((Budget-Beta1-c1*j)/c2) ;
        low_s2 = floor(max((Budget-Beta2-c1*j)/c2,1)) ;



        for k=low_s2:big_s2

            if k >= low_s2 && k<= up_s2

                tautotal = 0;
                totaltwo = 0;
                totalmobs = 0;
                average = 0;
                test = 0;
                samplese=0;


            for i=1:sample
                
                [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
                SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(j,k,[111112+i*4+t*100000 222223+i*4+t*100000 333334+i*4+t*100000 444444+i*4+t*100000],...
                warmup,m_0,...
                1e-3,Tol);

                if tau > 0
                    tautotal = tautotal+tau;
                    totaltwo = totaltwo+tau^2;

                    test = test+1;
                end

                totalmobs = totalmobs + mobs;
                
            end
            testnumber(j+(k-1)*big_s1)=test;
            if tautotal >0

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
    disp(['stage1 end mk=' num2str(m_0) 'totalm=' num2str(totalm)])

    sencond_tauavma=tauma;
    sencond_sema=sema;
    sencond_moma=moma;
    second_tuasqrma=tausqrma;
    second_test=testma;
    second_tautoma=tautoma;
    nose=0;
    Rinott = H_Rinotts(sample, 1-alpha/2, system_number);
    [best_solution]=NSGS(sencond_sema,sencond_tauavma,sample,Rinott);
    xlswrite('test11_testnumber',int2str(t),'.xls',testnumber);
end


