function [tau,AvgWaitingTime, Se_AvgWaitingTime, iseed, err,flag, ncalls_RF, Cost, SecurityLevel] = Bisection_Simulation(S1,S2,Epsilon,iseed0,...
    Num_Warmup,Num_Simulation,...
    Tol,...
    lb,ub)
    global log_test
    
    ncalls_RF = 0;
    CurUb = ub;
    CurLb = lb;
    
    [AvgWaitingTime,Se_AvgWaitingTime, iseed, R1, R2, P, Cost, SecurityLevel] = Simulation_AirportModel(CurUb,S1,S2,iseed0,...
                            Num_Warmup,Num_Simulation);

    ncalls_RF = ncalls_RF + Num_Simulation + Num_Warmup;
% 
%     AvgWT_lb = Simulation_AirportModel(CurLb,S1,S2,iseed0,...
%                         Num_Warmup,Num_Simulation);
% 
%     ncalls_search = ncalls_RF + Num_Simulation;
    
    while 1 == 1
        disp([ '(' num2str(CurLb) ',' num2str(CurUb) ')' ])
        fprintf(log_test,'\t\t\t (Lb, Ub) = (%f, %f)\n',CurLb, CurUb);
%         xdiff = CurUb - CurLb; % > 0
%         ydiff = AvgWaitingTime - AvgWT_lb; % < 0
%         guess = CurUb - (AvgWaitingTime-Epsilon/60)*(xdiff/ydiff);
        
        guess = (CurUb + CurLb)/2;

        [new_avgtime,new_se_AvgWaitingTime,iseed, R1, R2, P, Cost, SecurityLevel] = Simulation_AirportModel(guess,S1,S2,iseed0,...
                    Num_Warmup,Num_Simulation);

        ncalls_RF = ncalls_RF + Num_Simulation + Num_Warmup;

        if new_avgtime < Epsilon
            AvgWaitingTime = new_avgtime;
            Se_AvgWaitingTime = new_se_AvgWaitingTime;
            CurUb = guess;
        else
            CurLb = guess;
            AvgWT_lb = new_avgtime;
        end
        
        if CurUb - CurLb < Tol
            tau = CurUb;
            err = CurUb - CurLb;
            flag = 1;
            fprintf(log_test,'\t\t\t (Lb, Ub)(Terminated) = (%f, %f)\n',CurLb, CurUb);
            return
        end
    end
    
end