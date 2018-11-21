function [tau, meanWaitingTime, Se_AvgWaitingTime, flag, ncalls, iseed] = SolveTau_Simulation_Secant(Epsilon,S1,S2,Start,meanWaitingTime_tau_budget,iseedk,...
                Num_Warmup,Num_Simulation,...
                Delta0,Tol)
%  purpose: To solve the min possible tau which satisfies the mean waiting
%           time constraint when the min possible tau which satisfies the
%           budget does not satisfy the mean waiting time constraint.
%           The lower bound of the search region should be the min possible
%           tau which satisfies the budget but not the mean waiting time
%           constraint.
%           The upper bound of the search region should be a tau which
%           saisfies the mean waiting time constraint by guess or searching
%           using fminsearch().
%           After the lower bound and the upper bound are found, the
%           procedure apply the regular falsi method to solve the min
%           possible tau which satisfies both mean waiting time and budget
%           constraint.
%  input:
%   Epsilon:    The mean waiting time constraint.
%   S1:         The number of server in non-selectee lane
%   S2:         The number of server in selectee lane
%   Start:      The lower bound of the search region. This lower bound 
%               should be the min possible tau which satisfies 
%               the budget but not the mean waiting time constraint.
%  output:
%   tau:                The min possible tau which satisfies both mean waiting time
%                       and budget.
%   meanWaitingTime:    The corrsponded mean waiting time
%   flag:               The indicator of the validity of tau and
%                       meanWaitingTime. Valid if flag > 0 otherwise both
%                       the tau and meanWaitingTime are not valid.
    global log_test
    
    
    ncalls = 0;
    [ub, flag, meanWaitingTime_tau_guess, Se_AvgWaitingTime_guess] = FindUpperBound_Guess_Simulation(S1,S2,Epsilon,iseedk,Num_Warmup,Num_Simulation);
    ncalls = ncalls + Num_Simulation;
    
    NoInitialUb = 1;
    if flag < 0
        %guessed tau did not pass the WT criterion
%         if ub >= Start
%             %Start = ub; %update the lower bound
%         else
%             %Do nothing, keep the original lower bound
%         end
    else
        if ub < Start
            %irreasonable new tau (not satisfy the budget constraint)
        else
            NoInitialUb = 0;
            %Do nothing, use this upper bound directly
        end
    end
    
    if NoInitialUb == 1
        %Try to find ub anyway...
        fprintf(log_test,'\t\t\t Ub Initial Guess = %f (Infeasible)\n',ub);
        tau_guess = ub;
        tau_budget = Start;
        
        
        [lb, ub, flag, iseed, ncalls_search] = FindBoundingInterval(Epsilon,S1,S2,tau_budget,meanWaitingTime_tau_budget, tau_guess, meanWaitingTime_tau_guess,iseedk,...
                Num_Warmup,Num_Simulation,...
                Delta0,Tol,100);
        
        
%         [lb, ub, flag, meanWaitingTime, Se_AvgWaitingTime, iseed, ncalls_search] = FindUpperBound_Search_Simulation_Delta(S1,S2,Epsilon,iseedk,...
%                     Num_Warmup,Num_Simulation,...
%                     Start,1,Delta0); %Search the tau with AvgWT < Epsilon in (lb,1)
        ncalls = ncalls + ncalls_search;
        if flag < 0
            %No Solution...
            tau = [];
            meanWaitingTime = [];
            Se_AvgWaitingTime = [];
            flag = -1;
            return
        end
        Start = lb;
        fprintf(log_test,'\t\t\t Delta Approach (Lb, Ub) = (%f,%f)\n',Start, ub);
    else
%         guess = (Start + ub)/2;
%         [meanWaitingTime,Se_AvgWaitingTime] = Simulation_AirportModel(guess,S1,S2,iseedk,...
%                             Num_Warmup,Num_Simulation);
% 
%         ncalls = ncalls + Num_Simulation;
%        
%         
%         if meanWaitingTime*60 < Epsilon
%             ub = guess;
%         else
%             Start = guess;
%         end
%         [lb, ub_guess, flag, meanWaitingTime, Se_AvgWaitingTime, iseed, ncalls_search] = FindUpperBound_Search_Simulation_Delta(S1,S2,Epsilon,iseedk,...
%                     Num_Warmup,Num_Simulation,...
%                     Start,1,Delta0); %Search the tau with AvgWT < Epsilon in (lb,1)
%         ncalls = ncalls + ncalls_search;
%         
%         if lb > Start
%             Start = lb;
%         end
%         if ub_guess < ub
%             ub = ub_guess;
%         end
        fprintf(log_test,'\t\t\t Ub Initial Guess = %f (Feasible)\n',ub);
        fprintf(log_test,'\t\t\t Use Lb Directly (Lb, Ub) = (%f,%f)\n',Start, ub);
    end

                    
    [tau, meanWaitingTime, Se_AvgWaitingTime, iseed, err, flag, ncalls_RF] = RegularFalsi_SimpleBisection_Simulation(S1,S2,Epsilon,iseedk,...
                    Num_Warmup,Num_Simulation,...
                    Tol,...
                    Start, ub);
    ncalls = ncalls + ncalls_RF;
end