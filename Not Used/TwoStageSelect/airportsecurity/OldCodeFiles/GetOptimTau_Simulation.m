function [tau,MeanWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed] = GetOptimTau_Simulation(S1,S2,iseed0,...
    Num_Warmup,Num_Simulation,...
    Tol_Search,Tol_Solve,Tol)
%  purpose: To evaluate the min possible tau which satisfy both the budget
%  constraint and mean waiting time constraint
%  input
%   S1:             The number of servers in non-selectee lane.
%   S2:             The number of servers in selectee lane.
%  output
%   tau:            The min possible tau which satisfy the budget and mean
%                   waiting time constraint.
%   MeanWaitingTime:The corresponded mean waiting time of min possible tau.
%   MinTau_Budget:  The min possible tau which satisfy the budget
%                   constraint only.
%   flag:           The indicator of the validity of min possible tau.
%                   Valid if flag = 0. (-1 : infeasible)
%  variables:
%   Epsilon:        Mean waiting time constraint
%   RefMatrix:      The cache saves the min possible tau given server (S1,S2). 
%                   Note this matrix does not save the min possible tau
%                   given S1 = 0 or S2 = 0.

    global Epsilon
    global RefMatrix
    iseedk = iseed0(1:4);
    mobs = 0;
    tau = [];
    MeanWaitingTime = [];
    Se_AvgWaitingTime = [];
    MinTau_Budget = [];
    iseed = iseedk;
    flag2 = 1; %2: Done without solving; 1: Done, -1: Not Feasible (Epsilon), -2: Not Feasible (Budget)
    
    if S1 < 0 || S2 < 0
        flag2 = 1;
        return
    end
    
    if S1 == 0
        MinTau_Budget = GetMinTau(S1,S2);
        if ~isreal(MinTau_Budget) || MinTau_Budget > 0
            flag2 = 1;
            MinTau_Budget = [];
            return
        end

        [MeanWaitingTime,Se_AvgWaitingTime,iseed] = Simulation_AirportModel(0,S1,S2,iseedk,Num_Warmup,Num_Simulation);
        mobs = mobs + Num_Simulation;
        if MeanWaitingTime*60 < Epsilon
            tau = 0;
            flag2 = 0;
        else
            tau = [];
            flag2 = 1;
        end
        return
    end
    if S2 == 0
        MinTau_Budget = GetMinTau(S1,S2);
        if ~isreal(MinTau_Budget)
            flag2 = 1;
            MinTau_Budget = [];
            return
        end
        MinTau_Budget = 1;
        %iseed_vector = iseed0(1:4);
        [MeanWaitingTime,Se_AvgWaitingTime,iseed] = Simulation_AirportModel(1,S1,S2,iseedk,Num_Warmup,Num_Simulation);
        mobs = mobs + Num_Simulation;
        if MeanWaitingTime*60 < Epsilon
            tau = 1;
            flag2 = 0;
        else
            tau = [];
            flag2 = 1;
        end
        return
    end
    
    if S1 <= size(RefMatrix,1) && S2 <= size(RefMatrix,2)
        if RefMatrix(S1,S2) > 0 && S1 > 0 && S2 > 0
             %MinTau_Budget = GetMinTau(S1,S2);
             %tau = RefMatrix(S1,S2);
             %iseed_vector = iseed0(1:4);
             %MeanWaitingTime = Simulation_AirportModel(MinTau_Budget,S1,S2,iseed_vector,Num_Warmup,Num_Simulation)*60;
             %flag = 1;
             %return
        end
    end
    
    MinTau_Budget = GetMinTau(S1,S2);
    if ~isreal(MinTau_Budget)
        flag2 = 1;
        MinTau_Budget = [];
        return
    end
    if MinTau_Budget < 0
        MinTau_Budget = 0;
    end
    %iseed_vector = iseed0(1:(2+S1+S2));
    [MeanWaitingTime,Se_AvgWaitingTime,iseed] = Simulation_AirportModel(MinTau_Budget,S1,S2,iseedk,Num_Warmup,Num_Simulation);
    mobs = mobs + Num_Simulation;
    if MeanWaitingTime*60  < Epsilon
        tau = MinTau_Budget;
        if S1 > 0 && S2 > 0
            RefMatrix(S1,S2) = tau;
        end
        flag2 = 0;
        return
    else
        [tau , MeanWaitingTime, Se_AvgWaitingTime, flag, ncalls, iseed] = SolveTau_Simulation(Epsilon,S1,S2,MinTau_Budget,iseedk,...
            Num_Warmup,Num_Simulation,...
            Tol_Search,Tol_Solve,Tol);
        mobs = mobs + ncalls;
        if flag < 0
            flag2 = 1;
        else
            flag2 = 0;
        end
        if flag2 == 1
            tau = [];
            MeanWaitingTime = [];
        else
            if S1 > 0 && S2 > 0
                RefMatrix(S1,S2) = tau;
            end
        end
        
        return
    end
end
