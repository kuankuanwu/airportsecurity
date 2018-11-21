function [lb, ub, flag, meanWaitingTime, Se_AvgWaitingTime, iseed, ncalls_search] = FindUpperBound_Search_Simulation_Delta(S1,S2,Epsilon,iseed0,...
    Num_Warmup,Num_Simulation,...
    lb,ub,delta0)
    
    ncalls_search = 0;
    delta = delta0;
    initguess = lb;
    while 1 == 1
        new_guess_right = initguess + delta;
        if new_guess_right <= ub
            [new_avgtime,Se_AvgWaitingTime,iseed] = Simulation_AirportModel(new_guess_right,S1,S2,iseed0,...
                Num_Warmup,Num_Simulation);
            ncalls_search = ncalls_search + Num_Simulation;
            disp(['WT_R_Raw = ' num2str(new_avgtime)])
            disp(['WT_R = ' num2str(new_avgtime*60)])
            if new_avgtime*60 < Epsilon
                meanWaitingTime = new_avgtime;
                ub = new_guess_right;
                flag = 1;
                return
            else
                lb = new_guess_right;
            end
        else
            if new_guess_right == (initguess + delta)
                iseed = iseed0;
            end
            break;
        end
        delta = delta*2;
        
    end
    
    ub = [];
    flag = -1;
    meanWaitingTime = [];
    Se_AvgWaitingTime = [];
end