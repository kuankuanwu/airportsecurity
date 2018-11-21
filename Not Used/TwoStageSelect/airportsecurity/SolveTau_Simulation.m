function [tau, meanWaitingTime, Se_AvgWaitingTime, flag, ncalls, iseed] = SolveTau_Simulation(Epsilon,S1,S2,Start,iseedk,...
                Num_Warmup,Num_Simulation,...
                Tol_Search,Tol_Solve,Tol)
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
    ncalls = 0;
    [ub, flag, meanWaitingTime, Se_AvgWaitingTime] = FindUpperBound_Guess_Simulation(S1,S2,Epsilon,iseedk,Num_Warmup,Num_Simulation);
    ncalls = ncalls + Num_Simulation;
    if ub < Start
        flag = -1;
    end
    if flag < 0
        %Try to find ub anyway...
        [ub, flag, meanWaitingTime, Se_AvgWaitingTime, iseed, ncalls_search] = FindUpperBound_Search_Simulation(S1,S2,Epsilon,iseedk,...
                    Num_Warmup,Num_Simulation,...
                    Tol_Search,...
                    Start,1); %Search the tau with AvgWT < Epsilon in (lb,1)
        ncalls = ncalls + ncalls_search;
        if flag < 0
            %No Solution...
            tau = [];
            meanWaitingTime = [];
            Se_AvgWaitingTime = [];
            flag = -1;
            return
        end
    end
    [tau, meanWaitingTime, Se_AvgWaitingTime, iseed, err, flag, ncalls_RF] = RegularFalsi_Bisection_Simulation(S1,S2,Epsilon,iseedk,...
                    Num_Warmup,Num_Simulation,...
                    Tol_Solve, Tol,...
                    Start,ub);
    ncalls = ncalls + ncalls_RF;
end