
%% NSGS
function [best_sys]=NSGS(var_tau_matrix,tau_matrix,sample,system_number,alpha)
global big_s1
global big_s2
global c1
global c2
global Beta1
global Beta2
global Budget
Tol = 1e-5;
m_0 = 200;
warmup = 200;
t=tpdf(sample-1,(1-alpha/2)^(1/(system_number-1)));
%% find the best
best_tau=10000;
best_s=0;
for j=1:big_s1
        up_s2 =  floor((Budget-Beta1-c1*j)/c2) ;
        low_s2 = floor(max((Budget-Beta2-c1*j)/c2,1)) ;
        for k=low_s2:up_s2
            
            if tau_matrix(j+(k-1)*big_s1) < best_tau
                best_tau = tau_matrix (j+(k-1)*big_s1);
                best_s = j+(k-1)*big_s1;
            end
            
        end
end
W = 0;
same_system=[];
%% compare the system
for j=1:big_s1
    up_s2 =  floor((Budget-Beta1-c1*j)/c2) ;
    low_s2 = floor(max((Budget-Beta2-c1*j)/c2,1)) ;
    
    for k=low_s2:up_s2
        if    j+(k-1)*big_s1 == best_s
            continue
        end
        W = t * ((var_tau_matrix(j+(k-1)*big_s1)+var_tau_matrix(best_s))/sample)^0.5; 

        if -tau_matrix(j+(k-1)*big_s1) >= -tau_matrix(best_s) - W
           same_system =[same_system;[j,k,tau_matrix(j+(k-1)*big_s1)]];
           disp(['same system = (' num2str(j) ','  num2str(k) ')'])
        end
    end 
    
end
%% ranking

H = H_Rinotts(system_number, 1-alpha/2, sample);
    
[a,~]=size(same_system);
if a == 1
    best_sys = same_system;
    return
end    
for i=1:a
    number = same_system(i,1)+(same_system(i,2)-1)*big_s1;
    N = max(sample,ceil((H*(var_tau_matrix(number))^0.5/delta)^2));
    sum_tau = tau_matrix(number)*sample;
    
    for j=1:N-sample 
         [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
                SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(same_system(i,1),same_system(i,2),[111112+i*3+j*100000 222223+i*3+j*100000 333334+i*3+j*100000 444444+i*3+j*100000],...
                warmup,m_0,...
                1e-3,Tol);
        
        if tau > 0
          sum_tau = sum_tau + tau ;  
        end
        same_system(i,3)=sum_tau/N;
    end
end    

%% find the end system

best_sys = zeros(1,2);
best_tau = 100;

for i= 1:a
    
    if same_system(i,3) < best_tau
        best_tau = same_system(i,3);
        best_sys = same_system(i,1:2);
    end
end



end