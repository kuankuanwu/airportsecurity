function [tau, AvgWaitingTime, Se_AvgWaitingTime, flag, ncalls, iseed, Cost, SecurityLevel] = SolveTau_Simulation_Bisection(Epsilon,S1,S2,tau_budget,AvgWaitingTime_tau_budget,iseedk,...
                Num_Warmup,Num_Simulation,...
                Delta0,Tol)
%  purpose: To solve the min possible tau which satisfies the Avg waiting
%           time constraint when the min possible tau which satisfies the
%           budget does not satisfy the Avg waiting time constraint.
%           The lower bound of the search region should be the min possible
%           tau which satisfies the budget but not the Avg waiting time
%           constraint.
%           The upper bound of the search region should be a tau which
%           saisfies the Avg waiting time constraint by guess or searching
%           using fminsearch().
%           After the lower bound and the upper bound are found, the
%           procedure apply the regular falsi method to solve the min
%           possible tau which satisfies both Avg waiting time and budget
%           constraint.
%  input:
%   Epsilon:    The Avg waiting time constraint.
%   S1:         The number of server in non-selectee lane
%   S2:         The number of server in selectee lane
%   Start:      The lower bound of the search region. This lower bound 
%               should be the min possible tau which satisfies 
%               the budget but not the Avg waiting time constraint.
%  output:
%   tau:                The min possible tau which satisfies both Avg waiting time
%                       and budget.
%   AvgWaitingTime:    The corrsponded Avg waiting time
%   flag:               The indicator of the validity of tau and
%                       AvgWaitingTime. Valid if flag > 0 otherwise both
%                       the tau and AvgWaitingTime are not valid.
    global log_test
    global outter
    global k
    ncalls = 0;
    
  
    
    Cost = [];
    SecurityLevel = [];
    
    lb = tau_budget;
    lb_func = AvgWaitingTime_tau_budget;
    
    iseed = iseedk;
    
    [ub_guess, AvgWaitingTime_tau_guess, Se_AvgWaitingTime_guess] = UpperBoundGuess_Simulation(S1,S2,Epsilon,iseedk,Num_Warmup,Num_Simulation);
    ncalls = ncalls + Num_Simulation + Num_Warmup;
    
    
    if ub_guess > lb && AvgWaitingTime_tau_guess < Epsilon
        %Upper bound guess is feasible
        ub = ub_guess;
        ub_func = AvgWaitingTime_tau_guess;
        
        fprintf(log_test,'\t\t\t Ub Initial Guess = %f (Feasible)\n',ub);
        fprintf(log_test,'\t\t\t Use Lb Directly (Lb, Ub) = (%f,%f)\n',lb, ub);
    else
        %Upper bound guess is infeasible
        if AvgWaitingTime_tau_guess < Epsilon
            %ub_guess > lb && AvgWaitingTime_tau_guess*60 < Epsilon ==> No Solution...
            tau = [];
            AvgWaitingTime = [];
            Se_AvgWaitingTime = [];
            flag = -1;
            return
        end
        fprintf(log_test,'\t\t\t Ub Initial Guess = %f (Infeasible)\n',ub_guess);
        [tau_minWT,lb,ub, tau_minWT_func, lb_func, ub_func, iseed, ncalls_search] = FindMinAvgWait(Epsilon,S1,S2,tau_budget,AvgWaitingTime_tau_budget, ub_guess, AvgWaitingTime_tau_guess,iseedk,...
                Num_Warmup,Num_Simulation,...
                Tol,100*1.01^k); % here , change r to k
            
        ncalls = ncalls + ncalls_search;
        
        if tau_minWT_func > Epsilon || tau_minWT < lb
            %No Solution...
            tau = [];
            AvgWaitingTime = [];
            Se_AvgWaitingTime = [];
            flag = -1;
            return
        else
            ub = tau_minWT;
        end
        fprintf(log_test,'\t\t\t Min Avg WT Finding: (Lb, Ub) = (%f,%f)\n',lb, ub);
    end
                     
    [tau, AvgWaitingTime, Se_AvgWaitingTime, iseed, err, flag, ncalls_RF, Cost, SecurityLevel] = Bisection_Simulation(S1,S2,Epsilon,iseedk,...
                    Num_Warmup,Num_Simulation,...
                    Tol,...
                    1000*1.01^k,...   %here
                    lb, ub);
    ncalls = ncalls + ncalls_RF;
end