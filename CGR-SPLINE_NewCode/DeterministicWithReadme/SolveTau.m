function [tau, meanWaitingTime, flag] = SolveTau(Epsilon,S1,S2,Start)
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
    global chainflag
    global Tol_wait_fmin
    global Iter_fmin
    global Tol 
    options=optimset('MaxIter',Iter_fmin,'TolFun',Tol_wait_fmin,'TolX',Tol);
    [ub, flag] = FindUpperBound_Guess(S1,S2,Epsilon);
    if ub < Start
        if flag > 0
            tau = [];
            meanWaitingTime = [];
            flag = -1;
            return
        end
        flag = -1;
    end                                                                                                                                                                                                                                                         
    if flag < 0
%         [tau, meanWaitingTime] = fminsearch(@(tau)MeanWait_MatrixInverse(tau,S1,S2),ub);
        if chainflag==0
            [tau, meanWaitingTime] = fminsearch(@(tau)MeanWait_MatrixInverse(tau,S1,S2),ub,options);
        else
            [tau, meanWaitingTime] = fminsearch(@(tau)Deterministic_TestBlockTridiagonal_LargeMatrix2(tau,S1,S2),ub);
        end
        if meanWaitingTime > Epsilon || tau < Start
            tau = [];
            meanWaitingTime = [];
            flag = -1;
            return
        else
            ub = tau;
        end
    end
    [tau, meanWaitingTime, err, flag] = RegularFalsi_Bisection(Start,ub,S1,S2,Epsilon);
end