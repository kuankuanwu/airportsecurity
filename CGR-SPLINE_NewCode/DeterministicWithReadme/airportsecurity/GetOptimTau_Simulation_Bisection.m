function [tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel,Constraints] = GetOptimTau_Simulation_Bisection(S1,S2,iseed0,...
    Num_Warmup,Num_Simulation,...
    Delta0,Tol)
%  purpose: To evaluate the min possible tau which satisfy both the budget
%  constraint and Avg waiting time constraint
%  input
%   S1:             The number of servers in non-selectee lane.
%   S2:             The number of servers in selectee lane.
%  output
%   tau:            The min possible tau which satisfy the budget and Avg
%                   waiting time constraint.
%   AvgWaitingTime: The corresponded Avg waiting time of min possible tau.
%   MinTau_Budget:  The min possible tau which satisfy the budget
%                   constraint only.
%   flag:           The indicator of the validity of min possible tau.
%                   Valid if flag = 0. (-1 : infeasible)
%  variables:
%   Epsilon:        Avg waiting time constraint
%   RefMatrix:      The cache saves the min possible tau given server (S1,S2). 
%                   Note this matrix does not save the min possible tau
%                   given S1 = 0 or S2 = 0.
    global cal_tau
    global Epsilon
    global Budget
    global RefMatrix_tau
    global RefMatrix_AvgWT
    global RefMatrix_Se_AvgWT
    global RefMatrix_SecurityLevel
    global RefMatrix_Cost
    global RefMatrix_mobs
    global log_test
    iseedk = iseed0(1:4);
    mobs = 0;
    tau = [];
    AvgWaitingTime = [];
    Se_AvgWaitingTime = [];
    MinTau_Budget = [];
    SecurityLevel = [];
    Constraints = [];
    iseed = iseedk;
    flag2 = 1; %2: Done without solving; 1: Done, -1: Not Feasible (Epsilon), -2: Not Feasible (Budget)
    fprintf(log_test,'\t\tStart Solving Tau: ((S1,S2), mk, iseed) = ((%d, %d), %d, (%d, %d, %d))\n',S1, S2, Num_Simulation, iseed0(1), iseed0(2), iseed0(3));
    disp([ 'Start Solving Tau: ((S1,S2), mk, iseed) = ' num2str(S1) ',' num2str(S2) ','  num2str(Num_Simulation) ',' num2str(iseed0(1)) ',' num2str(iseed0(2)) ',' num2str(iseed0(3)) ])
    
    if S1 < 0 || S2 < 0
        flag2 = 1;
        return
    end
    
    if S1 == 0
        MinTau_Budget = GetMinTau(S1,S2);
        if ~isreal(MinTau_Budget) || MinTau_Budget > 0
            flag2 = 1;
            MinTau_Budget = [];
            fprintf(log_test,'\t\tEnd Solving Tau: INFEASIBLE\n');
            disp(['End Solving Tau: INFEASIBLE'])
            return
        end

        [AvgWaitingTime,Se_AvgWaitingTime,iseed,R1,R2,P,Cost,SecurityLevel] = Simulation_AirportModel(0,S1,S2,iseedk,Num_Warmup,Num_Simulation);
        mobs = mobs + Num_Simulation + Num_Warmup;
        if AvgWaitingTime < Epsilon
            tau = 0;
            flag2 = 0;
            fprintf(log_test,'\t\tEnd Solving Tau: ((S1,S2), tau, WT, SE_WT, ncalls) = ((%d, %d),%f,%f,%f, %d)\n',S1, S2, tau, AvgWaitingTime-Epsilon, Se_AvgWaitingTime, mobs);
            disp(['End Solving Tau: ((S1,S2), tau, WT, SE_WT, ncalls) = ' num2str(S1) ',' num2str(S2) ',' num2str(tau) ',' num2str(AvgWaitingTime-Epsilon) ','  num2str(Se_AvgWaitingTime) ',' num2str(mobs)])
            fprintf(log_test,'\t\tEnd Solving Tau: (SecurityLevel, Cost) = (%f, %f)\n', SecurityLevel, Cost-Budget);
            disp(['End Solving Tau: (SecurityLevel, Cost) = ' num2str(SecurityLevel) ','  num2str(Cost-Budget)])
            Constraints = [AvgWaitingTime-Epsilon Cost-Budget];
        else
            tau = [];
            flag2 = 1;
            fprintf(log_test,'\t\tEnd Solving Tau: INFEASIBLE\n');
            disp(['End Solving Tau: INFEASIBLE'])
        end
        return
    end
    if S2 == 0
        MinTau_Budget = GetMinTau(S1,S2);
        if ~isreal(MinTau_Budget)
            flag2 = 1;
            MinTau_Budget = [];
            fprintf(log_test,'\t\tEnd Solving Tau: INFEASIBLE\n');
            disp(['End Solving Tau: INFEASIBLE'])
            return
        end
        MinTau_Budget = 1;
        %iseed_vector = iseed0(1:4);
        [AvgWaitingTime,Se_AvgWaitingTime,iseed,R1,R2,P,Cost,SecurityLevel] = Simulation_AirportModel(1,S1,S2,iseedk,Num_Warmup,Num_Simulation);
        mobs = mobs + Num_Simulation + Num_Warmup;
        if AvgWaitingTime < Epsilon
            tau = 1;
            flag2 = 0;
            fprintf(log_test,'\t\tEnd Solving Tau: ((S1,S2), tau, WT, SE_WT, ncalls) = ((%d, %d),%f,%f,%f, %d)\n',S1, S2, tau, AvgWaitingTime-Epsilon, Se_AvgWaitingTime, mobs);
            disp(['End Solving Tau: ((S1,S2), tau, WT, SE_WT, ncalls) = ' num2str(S1) ',' num2str(S2) ',' num2str(tau) ',' num2str(AvgWaitingTime-Epsilon) ','  num2str(Se_AvgWaitingTime) ',' num2str(mobs)])
            fprintf(log_test,'\t\tEnd Solving Tau: (SecurityLevel, Cost) = (%f, %f)\n', SecurityLevel, Cost-Budget);
            disp(['End Solving Tau: (SecurityLevel, Cost) = ' num2str(SecurityLevel) ','  num2str(Cost-Budget)])
            Constraints = [AvgWaitingTime-Epsilon Cost-Budget];
        else
            tau = [];
            flag2 = 1;
            fprintf(log_test,'\t\tEnd Solving Tau: INFEASIBLE\n');
            disp(['End Solving Tau: INFEASIBLE'])
        end
        return
    end
    
    if S1 <= size(RefMatrix_tau,1) && S2 <= size(RefMatrix_tau,2)
        if RefMatrix_tau(S1,S2) ~= -1 && S1 > 0 && S2 > 0
             MinTau_Budget = GetMinTau(S1,S2);
             tau = RefMatrix_tau(S1,S2);
             AvgWaitingTime = RefMatrix_AvgWT(S1,S2);
             Se_AvgWaitingTime = RefMatrix_Se_AvgWT(S1,S2);
             SecurityLevel = RefMatrix_SecurityLevel(S1,S2);
             Cost = RefMatrix_Cost(S1,S2);
             Constraints = [AvgWaitingTime-Epsilon Cost-Budget];
             %mobs = RefMatrix_mobs(S1,S2);

             %iseed_vector = iseed0(1:4);
             %AvgWaitingTime = Simulation_AirportModel(MinTau_Budget,S1,S2,iseed_vector,Num_Warmup,Num_Simulation);
             flag2 = 0;
             fprintf(log_test,'\t\tEnd Solving Tau (Skip): ((S1,S2), tau, WT, SE_WT, ncalls) = ((%d, %d),%f,%f,%f, %d)\n',S1, S2, tau, AvgWaitingTime-Epsilon, Se_AvgWaitingTime, mobs);
             disp(['End Solving Tau(skip): ((S1,S2), tau, WT, SE_WT, ncalls) = ' num2str(S1) ',' num2str(S2) ',' num2str(tau) ',' num2str(AvgWaitingTime-Epsilon) ','  num2str(Se_AvgWaitingTime) ',' num2str(mobs)])
             fprintf(log_test,'\t\tEnd Solving Tau (Skip): (SecurityLevel, Cost) = (%f, %f)\n', SecurityLevel, Cost-Budget);
             disp(['End Solving Tau(skip): (SecurityLevel, Cost) = ' num2str(SecurityLevel) ','  num2str(Cost-Budget)])
             return
        end
    end
    
    MinTau_Budget = GetMinTau(S1,S2);
    if ~isreal(MinTau_Budget)
        flag2 = 1;
        MinTau_Budget = [];
        fprintf(log_test,'\t\tEnd Solving Tau: INFEASIBLE\n');
        disp(['End Solving Tau: INFEASIBLE'])
        return
    end
    if MinTau_Budget < 0
        MinTau_Budget = 0;
    end
    %iseed_vector = iseed0(1:(2+S1+S2));
    [AvgWaitingTime,Se_AvgWaitingTime,iseed,R1,R2,P,Cost,SecurityLevel] = Simulation_AirportModel(MinTau_Budget,S1,S2,iseedk,Num_Warmup,Num_Simulation);
    mobs = mobs + Num_Simulation + Num_Warmup;
    if AvgWaitingTime  < Epsilon
        tau = MinTau_Budget;
        
        cal_tau=[cal_tau;[S1,S2,tau,Num_Simulation]];
        
        if S1 > 0 && S2 > 0
            RefMatrix_tau(S1,S2) = tau;
            RefMatrix_AvgWT(S1,S2) = AvgWaitingTime;
            RefMatrix_Se_AvgWT(S1,S2) = Se_AvgWaitingTime;
            RefMatrix_SecurityLevel(S1,S2) = SecurityLevel;
            RefMatrix_Cost(S1,S2) = Cost;
            RefMatrix_mobs(S1,S2) = mobs;
        end
        flag2 = 0;
        fprintf(log_test,'\t\tEnd Solving Tau: ((S1,S2), tau, WT, SE_WT, ncalls) = ((%d, %d),%f,%f,%f, %d)\n',S1, S2, tau, AvgWaitingTime-Epsilon, Se_AvgWaitingTime, mobs);
        disp(['End Solving Tau: ((S1,S2), tau, WT, SE_WT, ncalls) = ' num2str(S1) ',' num2str(S2) ',' num2str(tau) ',' num2str(AvgWaitingTime-Epsilon) ','  num2str(Se_AvgWaitingTime) ',' num2str(mobs)])
        fprintf(log_test,'\t\tEnd Solving Tau: (SecurityLevel, Cost) = (%f, %f)\n', SecurityLevel, Cost-Budget);
        disp(['End Solving Tau: (SecurityLevel, Cost) = ' num2str(SecurityLevel) ','  num2str(Cost-Budget)])
        Constraints = [AvgWaitingTime-Epsilon Cost-Budget];
        return
    else
        [tau , AvgWaitingTime, Se_AvgWaitingTime, flag, ncalls, iseed, Cost, SecurityLevel] = SolveTau_Simulation_Bisection(Epsilon,S1,S2,MinTau_Budget,AvgWaitingTime,iseedk,...
            Num_Warmup,Num_Simulation,...
            Delta0,Tol);
        mobs = mobs + ncalls;
        if flag < 0
            flag2 = 1;
        else
            flag2 = 0;
        end
        if flag2 == 1
            tau = [];
            AvgWaitingTime = [];
            fprintf(log_test,'\t\tEnd Solving Tau: INFEASIBLE\n');
            disp(['End Solving Tau: INFEASIBLE'])
            return
        else
            if S1 > 0 && S2 > 0
                RefMatrix_tau(S1,S2) = tau;
                RefMatrix_AvgWT(S1,S2) = AvgWaitingTime;
                RefMatrix_Se_AvgWT(S1,S2) = Se_AvgWaitingTime;
                RefMatrix_SecurityLevel(S1,S2) = SecurityLevel;
                RefMatrix_Cost(S1,S2) = Cost;
                RefMatrix_mobs(S1,S2) = mobs;
            end
            fprintf(log_test,'\t\tEnd Solving Tau: ((S1,S2), tau, WT, SE_WT, ncalls) = ((%d, %d),%f,%f,%f, %d)\n',S1, S2, tau, AvgWaitingTime-Epsilon, Se_AvgWaitingTime, mobs);
            disp(['End Solving Tau: ((S1,S2), tau, WT, SE_WT, ncalls) = ' num2str(S1) ',' num2str(S2) ',' num2str(tau) ',' num2str(AvgWaitingTime-Epsilon) ','  num2str(Se_AvgWaitingTime) ',' num2str(mobs)])
            fprintf(log_test,'\t\tEnd Solving Tau: (SecurityLevel, Cost) = (%f, %f)\n', SecurityLevel, Cost-Budget);
            disp(['End Solving Tau: (SecurityLevel, Cost) = ' num2str(SecurityLevel) ','  num2str(Cost-Budget)])
            Constraints = [AvgWaitingTime-Epsilon Cost-Budget];
            return
        end
        
        
    end
end
