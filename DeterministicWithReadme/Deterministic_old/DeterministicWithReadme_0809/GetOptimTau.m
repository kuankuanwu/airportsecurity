function [tau,MeanWaitingTime,MinTau_Budget,flag] = GetOptimTau(S1,S2,ub)
%  purpose: To evaluate the min possible tau which satisfy both the budget
%  constraint and mean waiting time constraint
%  input
%   S1:             The number of servers in non-selectee lane.
%   S2:             The number of servers in selectee lane.
%   ub:             Upper bound of searching region. Usually is equal to 1.
%  output
%   tau:            The min possible tau which satisfy the budget and mean
%                   waiting time constraint.
%   MeanWaitingTime:The corresponded mean waiting time of min possible tau.
%   MinTau_Budget:  The min possible tau which satisfy the budget
%                   constraint only.
%   flag:           The indicator of the validity of min possible tau.
%                   Valid if flag > 0. 
%                   The details description of the value of flag:
%                   2: Done without solving, just lookup the cache which
%                   saves the results.
%                   1: Done and return the min possible tau successfully
%                   -1: Done, but the min possible tau does not satisfy the
%                   mean waiting time constraint.
%                   -2: Done, but the min possible tau does not satisfy the
%                   budget constraint.
%                   -3: Invalid S1 or S2
%  variables:
%   Epsilon:        Mean waiting time constraint
%   RefMatrix:      The cache saves the min possible tau given server (S1,S2). 
%                   Note this matrix does not save the min possible tau
%                   given S1 = 0 or S2 = 0.

    global Epsilon
    global RefMatrix
    global evaluation
    tau = [];
    MeanWaitingTime = [];
    MinTau_Budget = [];
    flag = -1; %2: Done without solving; 1: Done, -1: Not Feasible (Epsilon), -2: Not Feasible (Budget)
    
    if S1 < 0 || S2 < 0
        flag = -3;
        return
    end
    
    if S1 == 0
        MinTau_Budget = GetMinTau(S1,S2);
        if ~isreal(MinTau_Budget) || MinTau_Budget > 0
            flag = -2;
            MinTau_Budget = [];
            return
        end
        MeanWaitingTime = MeanWait_MatrixInverse(0,S1,S2)*60;
        if MeanWaitingTime < Epsilon
            tau = 0;
            flag = 1;
        else
            tau = [];
            flag = -1;
        end
        return
    end
    if S2 == 0
        MinTau_Budget = GetMinTau(S1,S2);
        if ~isreal(MinTau_Budget)
            flag = -2;
            MinTau_Budget = [];
            return
        end
        MinTau_Budget = 1;
        MeanWaitingTime = MeanWait_MatrixInverse(1,S1,S2)*60;
        if MeanWaitingTime < Epsilon
            tau = 1;
            flag = 1;
        else
            tau = [];
            flag = -1;
        end
        return
    end
    
    if  S1 > 0 && S2 > 0 && RefMatrix(S1,S2) ~= -1
        MinTau_Budget = GetMinTau(S1,S2);
        tau = RefMatrix(S1,S2);
        MeanWaitingTime = MeanWait_MatrixInverse(tau,S1,S2);
        flag = 1;
        return
    end
    
    MinTau_Budget = GetMinTau(S1,S2);
    if ~isreal(MinTau_Budget)
        flag = -2;
        MinTau_Budget = [];
        return
    end
    if MinTau_Budget < 0
        MinTau_Budget = 0;
    end
    MeanWaitingTime = MeanWait_MatrixInverse(MinTau_Budget,S1,S2)*60;
    if MeanWaitingTime  < Epsilon
        tau = MinTau_Budget;
        if S1 > 0 && S2 > 0
            RefMatrix(S1,S2) = tau;
            evaluation = evaluation+1;
        end
        flag = 2;
        return
    else
        [tau , MeanWaitingTime, flag] = SolveTau(Epsilon,S1,S2,MinTau_Budget);
        if flag < 0
            tau = [];
            MeanWaitingTime = [];
        else
            if S1 > 0 && S2 > 0
                RefMatrix(S1,S2) = tau;
                evaluation = evaluation+1;
            end
        end
        
        return
    end
end
