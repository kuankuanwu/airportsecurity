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
Budget = 95;'
Epsilon = 15;
global N
global lambda
global mu1
global mu2
global theta
global RefMatrix
global Lambda
N = 50;
lambda = 40/60;
Lambda = 40/60;
mu1 = 10/60;
mu2 = 1/60;
theta = 0.55;
global L1
global L2
global d1
global d2
d1 = 0.7;
d2 = 0.98;
global chainflag
chainflag=0;
global batch_num
batch_num = 20;
global b_meanq
global b_meanq2
b_meanq = [];
b_meanq2 = [];
Num_Warmup = 200;
 S1=4;
 S2=22;

 global Num_Simulation
 tau =0.281730242174354;
 iseed =  [1653381726,2882905302,1955172277,44444444];
 Num_Simulation = 54000 ;

 
  p = CalcP(tau,theta);
  ExpectedVisitors_NonSelecteeLane = p*Num_Simulation;
  ExpectedVisitors_SelecteeLane = (1-p)*Num_Simulation;
  if ExpectedVisitors_NonSelecteeLane >= 2
       m_eobm_NonSelecteeLane = floor(ExpectedVisitors_NonSelecteeLane/2);
  else
       m_eobm_NonSelecteeLane = -1;
  end
  if ExpectedVisitors_SelecteeLane >= 2
        m_eobm_SelecteeLane = floor(ExpectedVisitors_SelecteeLane/2);
  else
        m_eobm_SelecteeLane = -1;
  end
 
 
for i = 1:1
%     iseed = rand(1:4);
    [bsumq,all,selected,nonselected,MeanWaitingTimeSim,SE]=Simulation_AirportModel(tau,S1,S2,iseed,Num_Warmup,Num_Simulation);
    disp('obm')
    disp(obm(all,53981,length(all)));
    selected = obm(selected,length(selected)-m_eobm_SelecteeLane+1,length(selected));
    non = obm(nonselected,length(nonselected)- m_eobm_NonSelecteeLane+1,length(nonselected));
    disp('obm_p :');
    p_real = length(nonselected) / Num_Simulation ;
    disp((1-p_real)^2*selected + p_real^2*non);
    disp('nobm')
    disp(nobm(all,1000,length(all)));
end
    

function [VarAvgWait] = obm(sample,batch_size,total_size)
   

global b_meanq
   sum = 0;
   y_bar = mean(sample);
   batch_num = total_size - batch_size + 1;
   for i = 1:batch_num
       b_mean = mean(sample(i:i+batch_size-1));
       sum = sum + (b_mean-y_bar)^2;
       b_meanq = [b_meanq b_mean];
   end 
   VarAvgWait = batch_size/((total_size-batch_size+1)*(total_size-batch_size))*sum;
end

function [VarAvgWait] = nobm(sample,batch_size,total_size)
   
global batch_num
global b_meanq2
   sum = 0;
   y_bar = mean(sample);
   batch_num = ceil(total_size/batch_size);
   for i = 1:batch_num
       if i == batch_num
           b_mean = mean(sample((i-1)*batch_size +1 : total_size ));
       else
           b_mean = mean(sample((i-1)*batch_size +1 : i*batch_size ));
       end
       sum = sum + (b_mean-y_bar)^2;
       b_meanq2 = [b_meanq2 b_mean];
   end 
   VarAvgWait = batch_size/(batch_num-1)*sum;
end