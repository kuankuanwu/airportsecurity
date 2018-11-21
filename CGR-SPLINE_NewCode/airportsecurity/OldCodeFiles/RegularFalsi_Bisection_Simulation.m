function [tau,meanWaitingTime, Se_AvgWaitingTime, iseed, err,flag, ncalls_RF] = RegularFalsi_Bisection_Simulation(S1,S2,Epsilon,iseed0,...
    Num_Warmup,Num_Simulation,...
    Tol_Solve,Tol,...
    lb,ub)
    
    ncalls_RF = 0;
    CurUb = ub;
    while 1 == 1
        [ub, flag, meanWaitingTime, Se_AvgWaitingTime, iseed, ncalls_Search] = FindUpperBound_Search_Simulation(S1,S2,Epsilon,iseed0,...
                        Num_Warmup,Num_Simulation,...
                        Tol_Solve,...
                        lb,CurUb);
        ncalls_RF = ncalls_RF + ncalls_Search;
        if flag < 0
            %No solution better than given ub
            tau = CurUb;
            iseed_vector = iseed0(1:4);
            [meanWaitingTime, Se_AvgWaitingTime, iseed] = Simulation_AirportModel(tau,S1,S2,iseed_vector,...
                        Num_Warmup,Num_Simulation);
            err = ub - lb;
            flag = 2;
            return
        else
            %better tau exist
            CurUb = ub;
            if ub - lb < Tol
                %The searching interval < Tol ==> return ub as answer
                tau = ub;
                iseed_vector = iseed0(1:4);
                [meanWaitingTime, Se_AvgWaitingTime, iseed] = Simulation_AirportModel(tau,S1,S2,iseed_vector,...
                            Num_Warmup,Num_Simulation);
                err = ub - lb;
                flag = 1;
                return
            end
            disp([ '(lb,ub) = (' num2str(lb) ',' num2str(ub) ')'])
        end
    end
    
end