function [ub, flag, meanWaitingTime, Se_AvgWaitingTime, iseed, ncalls_search] = FindUpperBound_Search_Simulation(S1,S2,Epsilon,iseed0,...
    Num_Warmup,Num_Simulation,...
    Tol_Search,...
    lb,ub)
    
    ncalls_search = 0;
    
    Step = (ub - lb)/2;
    InitCutPoint = (ub + lb)/2; %2^0 points
    iseed_vector = iseed0(1:4);
    [MeanWaitingTime_InitCutPoint,Se_AvgWaitingTime_InitCutPoint,iseed] = Simulation_AirportModel(InitCutPoint,S1,S2,iseed_vector,...
        Num_Warmup,Num_Simulation);
    ncalls_search = ncalls_search + Num_Simulation;
    if MeanWaitingTime_InitCutPoint*60 < Epsilon
        ub = InitCutPoint;
        meanWaitingTime = MeanWaitingTime_InitCutPoint;
        Se_AvgWaitingTime = Se_AvgWaitingTime_InitCutPoint;
        flag = 1;
        return
    end
    PrevCutPoints = InitCutPoint;
    i = 1;
    while 1 == 1
        CutPoints = zeros(1,2^i);
        Step = Step / 2;
        if Step < Tol_Search
            break;
        end
        for j = 1 : length(PrevCutPoints)
            LeftCutPoint = PrevCutPoints(j) - Step;
            RightCutPoint = PrevCutPoints(j) + Step;
            CutPoints(2*(j-1)+1) = LeftCutPoint;
            CutPoints(2*(j-1)+2) = RightCutPoint;
            iseed_vector = iseed0(1:4);
            [MeanWaitingTime_LeftCutPoint,Se_AvgWaitingTime_LeftCutPoint,iseed] = Simulation_AirportModel(LeftCutPoint,S1,S2,iseed_vector,...
                Num_Warmup,Num_Simulation);
            ncalls_search = ncalls_search + Num_Simulation;
            disp(num2str(LeftCutPoint))
            if MeanWaitingTime_LeftCutPoint*60 < Epsilon
                ub = LeftCutPoint;
                meanWaitingTime = MeanWaitingTime_LeftCutPoint;
                Se_AvgWaitingTime = Se_AvgWaitingTime_LeftCutPoint;
                flag = 1;
                return
            end
            iseed_vector = iseed0(1:4);
            [MeanWaitingTime_RightCutPoint,Se_AvgWaitingTime_RightCutPoint,iseed] = Simulation_AirportModel(RightCutPoint,S1,S2,iseed_vector,...
                Num_Warmup,Num_Simulation);
            ncalls_search = ncalls_search + Num_Simulation;
            disp(num2str(RightCutPoint))
            if MeanWaitingTime_RightCutPoint*60 < Epsilon
                ub = RightCutPoint;
                meanWaitingTime = MeanWaitingTime_RightCutPoint;
                Se_AvgWaitingTime = Se_AvgWaitingTime_RightCutPoint;
                flag = 1;
                return
            end
        end
        PrevCutPoints = CutPoints;
        i = i + 1;
    end
    
    ub = [];
    flag = -1;
    meanWaitingTime = [];
    Se_AvgWaitingTime = [];
end