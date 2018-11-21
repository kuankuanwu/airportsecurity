function [ub, flag, meanWaitingTime, Se_AvgWaitingTime, iseed, ncalls_search] = FindNewUpperBound_Search_Simulation_Delta(S1,S2,Epsilon,iseed0,...
    Num_Warmup,Num_Simulation,...
    lb,ub,iter)

    ncalls_search = 0;
    meanWaitingTime = Simulation_AirportModel(ub,S1,S2,iseed0,...
                        Num_Warmup,Num_Simulation);

    ncalls_search = ncalls_search + Num_Simulation;
    
    AvgWT_lb = Simulation_AirportModel(lb,S1,S2,iseed0,...
                        Num_Warmup,Num_Simulation);

    ncalls_search = ncalls_search + Num_Simulation;
    
    xdiff = ub - lb; % > 0
    ydiff = meanWaitingTime - AvgWT_lb; % < 0
    guess = ub - (meanWaitingTime-Epsilon/60)*(xdiff/ydiff);
    
    [new_avgtime,Se_AvgWaitingTime,iseed] = Simulation_AirportModel(guess,S1,S2,iseed0,...
                Num_Warmup,Num_Simulation);
            
    ncalls_search = ncalls_search + Num_Simulation;
    
    if new_avgtime*60 < Epsilon
        meanWaitingTime = new_avgtime;
        ub = guess;
        flag = 1;
        return
    end
    
    PtVector = guess;
    NewPtVector = zeros(2,1);
    
    for i = 1 : iter
        for j = 1 : length(PtVector)
            if j == 1
                guess_left = (lb + PtVector(j))/2;
            else
                guess_left = (PtVector(j-1) + PtVector(j))/2;
            end
            
            [new_avgtime,Se_AvgWaitingTime,iseed] = Simulation_AirportModel(guess_left,S1,S2,iseed0,...
                    Num_Warmup,Num_Simulation);

            ncalls_search = ncalls_search + Num_Simulation;

            if new_avgtime*60 < Epsilon
                meanWaitingTime = new_avgtime;
                ub = guess_left;
                flag = 1;
                return
            end
            NewPtVector(2*(i-1)+1) = guess_left;
            
            if j == length(PtVector)
                guess_right = (PtVector(j) + ub)/2;
            else
                guess_right = (PtVector(j) + PtVector(j+1))/2;
            end
            [new_avgtime,Se_AvgWaitingTime,iseed] = Simulation_AirportModel(guess_right,S1,S2,iseed0,...
                    Num_Warmup,Num_Simulation);

            ncalls_search = ncalls_search + Num_Simulation;

            if new_avgtime*60 < Epsilon
                meanWaitingTime = new_avgtime;
                ub = guess_right;
                flag = 1;
                return
            end
            NewPtVector(2*(i-1)+2) = guess_right;
        end
        PtVector = NewPtVector;
    end
    
    ub = [];
    flag = -1;
    meanWaitingTime = [];
    Se_AvgWaitingTime = [];
end