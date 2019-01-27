function [waiting_time,SE,avg_waitingtime]
          = Simulation(tau,s1,s2,iseed,Num_Warmup,Num_Simulation)
          
disp(iseed)
global c1
global c2
global Beta1
global Beta2
global theta
global lambda
global mu1
global mu2
global N
    
