function [bsumq,waiting_time,...
          waiting_time_selected,...
          waiting_time_non_selected,...
          AvgWait,...
          SeAvgWait, ...
          iseed,...
          R1, ...
          R2, ...
          p_final, ...
          Cost, ...
          SecurityLevel] ...
          = Simulation_AirportModel(tau,s1,s2,iseed,Num_Warmup,Num_Simulation)
      
      %%非原版 不要用!!
    disp(iseed)
    %iseed: A vector contains random seeds
    %iseed(1): random seed for interarrival time
    %iseed(2): random seed for risk level
    %iseed(2 + (1~s1)): random seed for servers in non-selectee lane
    %iseed(end - s2 + 1 ~ end): random seed for servers in selectee
    %lane
    global c1
    global c2
    global Beta1
    global Beta2

    
    global theta
    global lambda
    global mu1
    global mu2
    global N
    
    
    global waiting_time
    global waiting_time_selected
    global waiting_time_non_selected
    waiting_time_non_selected = [];
    waiting_time_selected = [];
    waiting_time = [];
    
    
    global Sim_Time
    global Time_Last_Event
    global ServerStatus1 % 0==>idle, 1==>busy
    global ServerStatus2 % 0==>idle, 1==>busy
    global NextDepartureTime1
    global NextDepartureTime2
    global NextVisitTime
    global QueueNum_NonSelecteeLane
    global TimeArrival_NonSelecteeLane
    global ProcessTime_NonSelecteeLane
    global RiskValue_NonSelecteeLane
    %global TimeArrival_NonSelecteeLane_Processing
    global QueueNum_SelecteeLane
    global TimeArrival_SelecteeLane
    global ProcessTime_SelecteeLane
    global RiskValue_SelecteeLane
    global Test_ProcessTimeRecord
    Test_ProcessTimeRecord = [];
    %global TimeArrival_SelecteeLane_Processing
    
    %global Num_DepartedVisitors
    global Num_EndedWaitingVisitors
    global Num_CountedVisitors
    
    %global Num_DepartedVisitors_NonSelecteeLane
    global Num_EndedWaitingVisitors_NonSelecteeLane
    global Num_CountedVisitors_NonSelecteeLane
    
    %global Num_DepartedVisitors_SelecteeLane
    global Num_EndedWaitingVisitors_SelecteeLane
    global Num_CountedVisitors_SelecteeLane
    
    global m_eobm
    global queue_eobm
    global Sumy_WaitingTime
    global Sum_WaitingTime
    global Sum2_WaitingTime
    global Sumd_WaitingTime
    
    global m_eobm_NonSelecteeLane
    global queue_eobm_NonSelecteeLane
    global Sumy_WaitingTime_NonSelecteeLane
    global Sum_WaitingTime_NonSelecteeLane
    global Sum2_WaitingTime_NonSelecteeLane
    global Sumd_WaitingTime_NonSelecteeLane
    
    global m_eobm_SelecteeLane
    global queue_eobm_SelecteeLane
    global Sumy_WaitingTime_SelecteeLane
    global Sum_WaitingTime_SelecteeLane
    global Sum2_WaitingTime_SelecteeLane
    global Sumd_WaitingTime_SelecteeLane
    
    global Temp_Area_NonSelectee
    global Temp_Area_Selectee
    global Temp_SimDuration
    
    global RiskSum
    global RiskSum_NonSelecteeLane
    global RiskSum_SelecteeLane
    
    global Test_SelecteeWaitingTime
    Test_SelecteeWaitingTime = [];
    
    global Next_Event_Type %0==>EmptyQueue, 1==>Arrival, 2==>Departure
    
    global d1
    global d2
    
    
    global bsumq
    bsumq = [];
    m_eobm = max(Num_Simulation -19,1);

    p = CalcP(tau,theta);
    ExpectedVisitors_NonSelecteeLane = p*Num_Simulation;
    ExpectedVisitors_SelecteeLane = (1-p)*Num_Simulation;
    if ExpectedVisitors_NonSelecteeLane >= 2
        m_eobm_NonSelecteeLane = floor(ExpectedVisitors_NonSelecteeLane/2);
    else
        m_eobm_NonSelecteeLane = -1;
    end
    if ExpectedVisitors_SelecteeLane >= 2
        m_eobm_SelecteeLane = floor(ExpectedVisitors_SelecteeLane/2);
    else
        m_eobm_SelecteeLane = -1;
    end
    
    
    %Initialize
    iseed = Initialize(s1,s2,iseed);
    
    while(Num_EndedWaitingVisitors < Num_Warmup + Num_Simulation)
        
        timing
        
        update_time_avg_status
        
        if Next_Event_Type == 1
            iseed = arrive(tau,iseed,Num_Warmup);
            continue;
        elseif Next_Event_Type == 2
            iseed = depart(iseed,Num_Warmup);
            continue;
        end
    end   
    
    
    
    [ AvgWait,...
      VarAvgWait,...
      AvgWait_NonSelecteeLane,...
      VarAvgWait_NonSelecteeLane,...
      AvgWait_SelecteeLane,...
      VarAvgWait_SelecteeLane,...
      ] = obm;
  
    disp('obm_code :');
    disp(num2str(VarAvgWait))
    p_final = Num_CountedVisitors_NonSelecteeLane/Num_CountedVisitors;
    if (m_eobm_NonSelecteeLane > 0 && Num_CountedVisitors_NonSelecteeLane > 0) && ...
       (m_eobm_SelecteeLane > 0 && Num_CountedVisitors_SelecteeLane > 0)
        VarAvgWait = p_final^2*VarAvgWait_NonSelecteeLane + (1-p_final)^2*VarAvgWait_SelecteeLane;
    elseif m_eobm_NonSelecteeLane > 0 && Num_CountedVisitors_NonSelecteeLane > 0
        VarAvgWait = VarAvgWait_NonSelecteeLane;
    elseif m_eobm_SelecteeLane > 0 && Num_CountedVisitors_SelecteeLane > 0
        VarAvgWait = VarAvgWait_SelecteeLane;
    end
    
    disp('obm_p_code :');
    disp(num2str(VarAvgWait))
    
    R1 = RiskSum_NonSelecteeLane/RiskSum;
    R2 = RiskSum_SelecteeLane/RiskSum;
        
    SecurityLevel = R1*d1 + R2*d2;
    Cost = c1*s1 + c2*s2 + CalcP(tau,theta)*Beta1 + (1-CalcP(tau,theta))*Beta2;
    AvgWait = (AvgWait_NonSelecteeLane + 1/mu1)*p_final + (AvgWait_SelecteeLane + 1/mu2)*(1-p_final);
%     waiting time ing queue 
    AvgWait = (AvgWait_NonSelecteeLane )*p_final + (AvgWait_SelecteeLane )*(1-p_final);
    SeAvgWait = sqrt(VarAvgWait);
   

    
end

function iseed = arrive(tau,iseed,Num_Warmup)
    global theta
    global lambda
    
    global Sim_Time
	global NextVisitTime
    % Schedule next arrival.
    [Interarrival,iseed(1)] = rexpo(1/lambda,iseed(1));
    NextVisitTime = Sim_Time + Interarrival;
    [RiskValue, iseed(2)] = Risk_Visitor(theta,iseed(2));
    
    if RiskValue <= tau
        iseed = arrive_non_selectee(iseed,Num_Warmup,RiskValue);
    else
        iseed = arrive_selectee(iseed,Num_Warmup,RiskValue);
    end
    
end

function iseed = arrive_non_selectee(iseed,Num_Warmup,RiskValue)
    global mu1
    global mu2
    global N
    global waiting_time_non_selected
    global waiting_time
    global Sim_Time
    global ServerStatus1 % 0==>idle, 1==>busy
    global NextDepartureTime1
    global QueueNum_NonSelecteeLane
    global QueueNum_SelecteeLane
    global TimeArrival_NonSelecteeLane
    global ProcessTime_NonSelecteeLane
    global RiskValue_NonSelecteeLane
    global Test_ProcessTimeRecord
    %global TimeArrival_NonSelecteeLane_Processing
    
    global Num_EndedWaitingVisitors
    global Num_CountedVisitors
    
    global Num_EndedWaitingVisitors_NonSelecteeLane
    global Num_CountedVisitors_NonSelecteeLane
    
    IdleServerIdx = 0;
    for i = 1 : length(ServerStatus1)
        if ServerStatus1(i) == 0
            IdleServerIdx = i;
        end
    end
    [ProcessTime,iseed(3)] = rexpo(1/mu1,iseed(3));
    %[Dummy,iseed(4)] = rexpo(1/mu2,iseed(4));
    Test_ProcessTimeRecord = [Test_ProcessTimeRecord ProcessTime];
    if IdleServerIdx == 0
        %Server Status are busy
        if QueueNum_NonSelecteeLane + QueueNum_SelecteeLane < N
			%if the sum of two queue lengths is less than N, add this visitor into queue.
            QueueNum_NonSelecteeLane = QueueNum_NonSelecteeLane + 1;
            TimeArrival_NonSelecteeLane = [TimeArrival_NonSelecteeLane Sim_Time];
            ProcessTime_NonSelecteeLane = [ProcessTime_NonSelecteeLane ProcessTime];
            RiskValue_NonSelecteeLane =   [RiskValue_NonSelecteeLane RiskValue];
        else
            %Visitor Balked (Overflowed)
        end
    else
        %Server is idle, let visitor go to server directly
        %Record the delay time
        DelayTime = 0;
        if Num_EndedWaitingVisitors >= Num_Warmup
            Num_CountedVisitors = Num_CountedVisitors + 1;
            Num_CountedVisitors_NonSelecteeLane = Num_CountedVisitors_NonSelecteeLane + 1;
            update_sum_variables(DelayTime,RiskValue)
            update_sum_variables_NonSelecteeLane(DelayTime,RiskValue)
            waiting_time_non_selected = [waiting_time_non_selected  DelayTime];
            waiting_time = [waiting_time  DelayTime];
        end
        Num_EndedWaitingVisitors = Num_EndedWaitingVisitors + 1;
        Num_EndedWaitingVisitors_NonSelecteeLane = Num_EndedWaitingVisitors_NonSelecteeLane + 1;
        
        ServerStatus1(IdleServerIdx) = 1;
        
        %ProcessTime = 1/mu1;
        NextDepartureTime1(IdleServerIdx) = Sim_Time + ProcessTime;
    end
end

function iseed = arrive_selectee(iseed,Num_Warmup,RiskValue)
    global mu2
    global mu1
    global N
    
    global waiting_time
    global waiting_time_selected
    global Sim_Time
    global ServerStatus2 % 0==>idle, 1==>busy
    global NextDepartureTime2
    global QueueNum_NonSelecteeLane
    global QueueNum_SelecteeLane
    global TimeArrival_SelecteeLane
    global ProcessTime_SelecteeLane
    global RiskValue_SelecteeLane
    global Test_ProcessTimeRecord
    %global TimeArrival_SelecteeLane_Processing
    
    global Num_EndedWaitingVisitors
    global Num_CountedVisitors
    
    global Num_EndedWaitingVisitors_SelecteeLane
    global Num_CountedVisitors_SelecteeLane
    
    [ProcessTime,iseed(3)] = rexpo(1/mu2,iseed(3));
    %[Dummy,iseed(3)] = rexpo(1/mu1,iseed(3));
    Test_ProcessTimeRecord = [Test_ProcessTimeRecord ProcessTime];
    if QueueNum_SelecteeLane == 0
        %The queue is empty
        
        IdleServerIdx = 0; %To find the ID of the first idle server (0 means all servers are busy)
        for i = 1 : length(ServerStatus2)
            if ServerStatus2(i) == 0
                IdleServerIdx = i; % The ith server is idle
                break;
            end
        end
        if IdleServerIdx == 0
            %Server Status are busy
            if QueueNum_NonSelecteeLane + QueueNum_SelecteeLane < N
                QueueNum_SelecteeLane = QueueNum_SelecteeLane + 1;
                TimeArrival_SelecteeLane = [TimeArrival_SelecteeLane Sim_Time];
                ProcessTime_SelecteeLane = [ProcessTime_SelecteeLane ProcessTime];
                RiskValue_SelecteeLane = [RiskValue_SelecteeLane RiskValue];
            else
                %Visitor Balked (Overflowed)
            end
        else
            %Server is idle, let visitor go to server directly
            %Record the delay time
            DelayTime = 0;
            if Num_EndedWaitingVisitors >= Num_Warmup
                Num_CountedVisitors = Num_CountedVisitors + 1;
                Num_CountedVisitors_SelecteeLane = Num_CountedVisitors_SelecteeLane + 1;
                update_sum_variables(DelayTime,RiskValue)
                update_sum_variables_SelecteeLane(DelayTime,RiskValue)
                waiting_time_selected = [waiting_time_selected ,DelayTime];
                waiting_time = [waiting_time ,DelayTime];
                
            end
            Num_EndedWaitingVisitors = Num_EndedWaitingVisitors + 1;
            Num_EndedWaitingVisitors_SelecteeLane = Num_EndedWaitingVisitors_SelecteeLane + 1;
            
            ServerStatus2(IdleServerIdx) = 1;
            
            %ProcessTime = 1/mu2;
            NextDepartureTime2(IdleServerIdx) = Sim_Time + ProcessTime;
        end
    else
        if QueueNum_NonSelecteeLane + QueueNum_SelecteeLane < N
            QueueNum_SelecteeLane = QueueNum_SelecteeLane + 1;
            TimeArrival_SelecteeLane = [TimeArrival_SelecteeLane Sim_Time];
            ProcessTime_SelecteeLane = [ProcessTime_SelecteeLane ProcessTime];
            RiskValue_SelecteeLane = [RiskValue_SelecteeLane RiskValue];
        else
            %Visitor Balked (Overflowed)
        end
    end
end


function iseed = depart(iseed,Num_Warmup)
    global NextDepartureTime1
    global NextDepartureTime2
    Next_Departure_Type = 0; %1==>NonSelectee %2==>Selectee
    Next_Departure_ServerIdx = 0;
    Min_Time_Next_Departure = 1.0e29;
    %Check next departure (Non selectee)
    if ~isempty(NextDepartureTime1)
        for j = 1 : length(NextDepartureTime1)
            if NextDepartureTime1(j) >= 0
                if NextDepartureTime1(j) < Min_Time_Next_Departure
                    Next_Departure_Type = 1;
                    Next_Departure_ServerIdx = j;
                    Min_Time_Next_Departure = NextDepartureTime1(j);
                end
            end
        end
    end
    %Check next departure (Selectee)
    if ~isempty(NextDepartureTime2)
        for j = 1 : length(NextDepartureTime2)
            if NextDepartureTime2(j) >= 0
                if NextDepartureTime2(j) < Min_Time_Next_Departure
                    Next_Departure_Type = 2;
                    Next_Departure_ServerIdx = j;
                    Min_Time_Next_Departure = NextDepartureTime2(j);
                end
            end
        end
    end
    
    if Next_Departure_Type == 1
        iseed = depart_non_selectee(iseed,Next_Departure_ServerIdx,Num_Warmup);
    else
        iseed = depart_selectee(iseed,Next_Departure_ServerIdx,Num_Warmup);
    end
end

function iseed = depart_non_selectee(iseed,Next_Departure_ServerIdx,Num_Warmup)
    global mu1
    global waiting_time_non_se
    global Sim_Time
    global ServerStatus1
    global NextDepartureTime1
    global QueueNum_NonSelecteeLane
    global TimeArrival_NonSelecteeLane
    global ProcessTime_NonSelecteeLane
    global RiskValue_NonSelecteeLane
    %global TimeArrival_NonSelecteeLane_Processing
    %global Num_DepartedVisitors
    global Num_EndedWaitingVisitors
    global Num_CountedVisitors
    %global Num_DepartedVisitors_NonSelecteeLane
    global Num_EndedWaitingVisitors_NonSelecteeLane
    global Num_CountedVisitors_NonSelecteeLane
    
    global waiting_time
   
    global waiting_time_selected
    
    global waiting_time_non_selected

    
    %Process the next visitor if the queue is not empty:
    if QueueNum_NonSelecteeLane == 0
        %Server status ==> idle
        ServerStatus1(Next_Departure_ServerIdx) = 0;
        NextDepartureTime1(Next_Departure_ServerIdx) = -1;
    else
        %Process the next waiting visitors
        
        %Record the delay time
        DelayTime = Sim_Time - TimeArrival_NonSelecteeLane(1);
        RiskValue = RiskValue_NonSelecteeLane(1);
        if Num_EndedWaitingVisitors >= Num_Warmup
            Num_CountedVisitors = Num_CountedVisitors + 1;
            Num_CountedVisitors_NonSelecteeLane = Num_CountedVisitors_NonSelecteeLane + 1;
            update_sum_variables(DelayTime,RiskValue)
            update_sum_variables_NonSelecteeLane(DelayTime,RiskValue)
            waiting_time_non_selected = [waiting_time_non_selected  DelayTime];
            waiting_time = [waiting_time  DelayTime];
        end

        Num_EndedWaitingVisitors = Num_EndedWaitingVisitors + 1;
        Num_EndedWaitingVisitors_NonSelecteeLane = Num_EndedWaitingVisitors_NonSelecteeLane + 1;
        
        QueueNum_NonSelecteeLane = QueueNum_NonSelecteeLane - 1;
        TimeArrival_NonSelecteeLane(1) = [];
        ProcessTime = ProcessTime_NonSelecteeLane(1);
        ProcessTime_NonSelecteeLane(1) = [];
        RiskValue_NonSelecteeLane(1) = [];
        %ProcessTime = 1/mu1;
        NextDepartureTime1(Next_Departure_ServerIdx) = Sim_Time + ProcessTime;
        
        
    end
end

function iseed = depart_selectee(iseed,Next_Departure_ServerIdx,Num_Warmup)
    global mu2
    
    global Sim_Time
    global ServerStatus2
    global NextDepartureTime2
    global QueueNum_SelecteeLane
    global TimeArrival_SelecteeLane
    global ProcessTime_SelecteeLane
    global RiskValue_SelecteeLane
    %global TimeArrival_SelecteeLane_Processing
    %global Num_DepartedVisitors
    global Num_EndedWaitingVisitors
    global Num_CountedVisitors
    %global Num_DepartedVisitors_SelecteeLane
    global Num_EndedWaitingVisitors_SelecteeLane
    global Num_CountedVisitors_SelecteeLane
    
    global waiting_time
   
    global waiting_time_selected
    
    global waiting_time_non_selected
    
    %Process the next visitor if the queue is not empty:
    if QueueNum_SelecteeLane == 0
        %Server status ==> idle
        ServerStatus2(Next_Departure_ServerIdx) = 0;
        NextDepartureTime2(Next_Departure_ServerIdx) = -1;
    else
        %Process the next waiting visitors
        %Record the delay time
        DelayTime = Sim_Time - TimeArrival_SelecteeLane(1);
        RiskValue = RiskValue_SelecteeLane(1);
        if Num_EndedWaitingVisitors >= Num_Warmup
            Num_CountedVisitors = Num_CountedVisitors + 1;
            Num_CountedVisitors_SelecteeLane = Num_CountedVisitors_SelecteeLane + 1;
            update_sum_variables(DelayTime,RiskValue)
            update_sum_variables_SelecteeLane(DelayTime,RiskValue)
            waiting_time_selected = [waiting_time_selected  DelayTime];
            waiting_time = [waiting_time  DelayTime];
        end
        Num_EndedWaitingVisitors = Num_EndedWaitingVisitors + 1;
        Num_EndedWaitingVisitors_SelecteeLane = Num_EndedWaitingVisitors_SelecteeLane + 1;
    
        QueueNum_SelecteeLane = QueueNum_SelecteeLane - 1;
        TimeArrival_SelecteeLane(1) = [];
        ProcessTime = ProcessTime_SelecteeLane(1);
        ProcessTime_SelecteeLane(1) = [];
        RiskValue_SelecteeLane(1) = [];
        %ProcessTime = 1/mu2;
        NextDepartureTime2(Next_Departure_ServerIdx) = Sim_Time + ProcessTime;
    end
end


function iseed = Initialize(s1,s2,iseed)
    global lambda
    
    global Sim_Time
    global Time_Last_Event
    global ServerStatus1
    global ServerStatus2
    global NextDepartureTime1
    global NextDepartureTime2
    global NextVisitTime
    global QueueNum_NonSelecteeLane
    global TimeArrival_NonSelecteeLane
    global ProcessTime_NonSelecteeLane
    global RiskValue_NonSelecteeLane
    global TimeArrival_NonSelecteeLane_Processing
    global QueueNum_SelecteeLane
    global TimeArrival_SelecteeLane
    global ProcessTime_SelecteeLane
    global RiskValue_SelecteeLane
    global TimeArrival_SelecteeLane_Processing
    
    
    
    %global Num_DepartedVisitors
    global Num_EndedWaitingVisitors
    global Num_CountedVisitors
    
    %global Num_DepartedVisitors_NonSelecteeLane
    global Num_EndedWaitingVisitors_NonSelecteeLane
    global Num_CountedVisitors_NonSelecteeLane
    
    %global Num_DepartedVisitors_SelecteeLane
    global Num_EndedWaitingVisitors_SelecteeLane
    global Num_CountedVisitors_SelecteeLane
    
    global queue_eobm
    global queue_eobm_NonSelecteeLane
    global queue_eobm_SelecteeLane
    global Sumy_WaitingTime
    global Sumy_WaitingTime_NonSelecteeLane
    global Sumy_WaitingTime_SelecteeLane
    global Sum_WaitingTime
    global Sum_WaitingTime_NonSelecteeLane
    global Sum_WaitingTime_SelecteeLane
    global Sum2_WaitingTime
    global Sum2_WaitingTime_NonSelecteeLane
    global Sum2_WaitingTime_SelecteeLane
    global Sumd_WaitingTime
    global Sumd_WaitingTime_NonSelecteeLane
    global Sumd_WaitingTime_SelecteeLane
    
    global Temp_Area_NonSelectee
    global Temp_Area_Selectee
    global Temp_SimDuration
    
    global RiskSum
    global RiskSum_NonSelecteeLane
    global RiskSum_SelecteeLane
    
    
    Sim_Time = 0;
    Time_Last_Event = 0;
    ServerStatus1 = zeros(1,s1); % 0==>idle, 1==>busy
    ServerStatus2 = zeros(1,s2); % 0==>idle, 1==>busy
    NextDepartureTime1 = -ones(1,s1);
    NextDepartureTime2 = -ones(1,s2);
    QueueNum_NonSelecteeLane = 0;
    TimeArrival_NonSelecteeLane = [];
    ProcessTime_NonSelecteeLane = [];
    RiskValue_NonSelecteeLane = [];
    TimeArrival_NonSelecteeLane_Processing = -ones(1,s1);
    QueueNum_SelecteeLane = 0;
    TimeArrival_SelecteeLane = [];
    ProcessTime_SelecteeLane = [];
    RiskValue_SelecteeLane = [];
    TimeArrival_SelecteeLane_Processing = -ones(1,s2);
    
    
    Num_DepartedVisitors = 0;
    Num_EndedWaitingVisitors = 0;
    Num_CountedVisitors = 0;
    Num_DepartedVisitors_NonSelecteeLane = 0;
    Num_EndedWaitingVisitors_NonSelecteeLane = 0;
    Num_CountedVisitors_NonSelecteeLane = 0;
    Num_DepartedVisitors_SelecteeLane = 0;
    Num_EndedWaitingVisitors_SelecteeLane = 0;
    Num_CountedVisitors_SelecteeLane = 0;
    
    Sum_WaitingTime = 0;
    Sum_WaitingTime_NonSelecteeLane = 0;
    Sum_WaitingTime_SelecteeLane = 0;
    queue_eobm = [];
    queue_eobm_NonSelecteeLane = [];
    queue_eobm_SelecteeLane = [];
    Sumy_WaitingTime = 0;
    Sumy_WaitingTime_NonSelecteeLane = 0;
    Sumy_WaitingTime_SelecteeLane = 0;
    Sum2_WaitingTime = 0;
    Sum2_WaitingTime_NonSelecteeLane = 0;
    Sum2_WaitingTime_SelecteeLane = 0;
    Sumd_WaitingTime = 0;
    Sumd_WaitingTime_NonSelecteeLane = 0;
    Sumd_WaitingTime_SelecteeLane = 0;
    
    Temp_Area_NonSelectee = 0;
    Temp_Area_Selectee = 0;
    Temp_SimDuration = 0;
    
    RiskSum = 0;
    RiskSum_NonSelecteeLane = 0;
    RiskSum_SelecteeLane = 0;
    
   [Interarrival,iseed(1)] = rexpo(1/lambda,iseed(1));
   NextVisitTime = Sim_Time + Interarrival;
   
end

function timing
    global Sim_Time
    global NextDepartureTime1
    global NextDepartureTime2
    global NextVisitTime
    
    global Next_Event_Type %0==>EmptyQueue, 1==>Arrival, 2==>Departure
    Next_Event_Type = 0;
    Min_Time_Next_Event = 1.0e29;
	%Determine the event type of the next event to occur.
    %Check next arrival
    if NextVisitTime < Min_Time_Next_Event
        Next_Event_Type = 1;
        Min_Time_Next_Event = NextVisitTime;
    end
    %Check next departure (Non selectee)
    if ~isempty(NextDepartureTime1)
        for j = 1 : length(NextDepartureTime1)
            if NextDepartureTime1(j) >= 0
                if NextDepartureTime1(j) < Min_Time_Next_Event
                    Next_Event_Type = 2;
                    Min_Time_Next_Event = NextDepartureTime1(j);
                end
            end
        end
    end
    %Check next departure (Selectee)
    if ~isempty(NextDepartureTime2)
        for j = 1 : length(NextDepartureTime2)
            if NextDepartureTime2(j) >= 0
                if NextDepartureTime2(j) < Min_Time_Next_Event
                    Next_Event_Type = 2;
                    Min_Time_Next_Event = NextDepartureTime2(j);
                end
            end
        end
    end
    
    Sim_Time = Min_Time_Next_Event;
end

function update_time_avg_status
    global Sim_Time
    
    global Time_Last_Event
    
    global ServerStatus1 % 0==>idle, 1==>busy
    global ServerStatus2 % 0==>idle, 1==>busy
    
    global QueueNum_NonSelecteeLane
    global QueueNum_SelecteeLane
    
    global Temp_Area_NonSelectee
    global Temp_Area_Selectee
    
    global Num_DepartedVisitors
    global Temp_SimDuration
    
    Time_Since_Last_Event = Sim_Time - Time_Last_Event;
    if Num_DepartedVisitors > 200
        Temp_SimDuration = Temp_SimDuration + Time_Since_Last_Event;
        Temp_Area_NonSelectee = Temp_Area_NonSelectee + (QueueNum_NonSelecteeLane + sum(ServerStatus1))*Time_Since_Last_Event;
        Temp_Area_Selectee = Temp_Area_Selectee + (QueueNum_SelecteeLane + sum(ServerStatus2))*Time_Since_Last_Event;
    end
    Time_Last_Event = Sim_Time;
end

function [alpha,iseed] = Risk_Visitor(theta,iseed)
    %Truncated exponential Distribution: Get the fractional part of
    %exponential random variable
    [Tmp,iseed] = rexpo(theta,iseed);
    Tmp = Tmp - floor(Tmp);
    alpha = Tmp;
    return 
end

function update_sum_variables(WaitingTime,RiskValue)
    
    global m_eobm
    global queue_eobm
    global Sumy_WaitingTime
    global Sum_WaitingTime
    global Sum2_WaitingTime
    global Sumd_WaitingTime
    global Num_CountedVisitors
    
    global RiskSum
    global bsumq
    
    RiskSum = RiskSum + RiskValue;
    
    if Num_CountedVisitors <= m_eobm
        %process the first m observations---the first complete batch     
        Sumy_WaitingTime = Sumy_WaitingTime + WaitingTime;
         queue_eobm = [queue_eobm WaitingTime];
        if Num_CountedVisitors == m_eobm
            Sum_WaitingTime = Sumy_WaitingTime;
            Sum2_WaitingTime = Sum_WaitingTime*Sum_WaitingTime;
            Sumd_WaitingTime = 0;
        end
    else
        %process observations m+1 through n
        Sumd_WaitingTime = Sumd_WaitingTime + queue_eobm(1);
        Sumy_WaitingTime = Sumy_WaitingTime + WaitingTime;
        bsum = Sumy_WaitingTime - Sumd_WaitingTime;
        Sum_WaitingTime = Sum_WaitingTime + bsum;
        Sum2_WaitingTime = Sum2_WaitingTime + bsum*bsum;
        queue_eobm(1) = [];
        queue_eobm = [queue_eobm WaitingTime];
        bsumq = [bsumq bsum];
    end
end

function update_sum_variables_NonSelecteeLane(WaitingTime,RiskValue)
    
    global m_eobm_NonSelecteeLane
    global queue_eobm_NonSelecteeLane
    global Sumy_WaitingTime_NonSelecteeLane
    global Sum_WaitingTime_NonSelecteeLane
    global Sum2_WaitingTime_NonSelecteeLane
    global Sumd_WaitingTime_NonSelecteeLane
    global Num_CountedVisitors_NonSelecteeLane
    
    global RiskSum_NonSelecteeLane
    
    
    
    RiskSum_NonSelecteeLane = RiskSum_NonSelecteeLane + RiskValue;
    if m_eobm_NonSelecteeLane < 0
        return
    end
    
    if Num_CountedVisitors_NonSelecteeLane <= m_eobm_NonSelecteeLane
        %process the first m observations---the first complete batch     
        Sumy_WaitingTime_NonSelecteeLane = Sumy_WaitingTime_NonSelecteeLane + WaitingTime;
        queue_eobm_NonSelecteeLane = [queue_eobm_NonSelecteeLane WaitingTime];
        if Num_CountedVisitors_NonSelecteeLane == m_eobm_NonSelecteeLane
            Sum_WaitingTime_NonSelecteeLane = Sumy_WaitingTime_NonSelecteeLane;
            Sum2_WaitingTime_NonSelecteeLane = Sum_WaitingTime_NonSelecteeLane*Sum_WaitingTime_NonSelecteeLane;
            Sumd_WaitingTime_NonSelecteeLane = 0;
        end
    else
        %process observations m+1 through n
        Sumd_WaitingTime_NonSelecteeLane = Sumd_WaitingTime_NonSelecteeLane + queue_eobm_NonSelecteeLane(1);
        Sumy_WaitingTime_NonSelecteeLane = Sumy_WaitingTime_NonSelecteeLane + WaitingTime;
        bsum = Sumy_WaitingTime_NonSelecteeLane - Sumd_WaitingTime_NonSelecteeLane;
        Sum_WaitingTime_NonSelecteeLane = Sum_WaitingTime_NonSelecteeLane + bsum;
        Sum2_WaitingTime_NonSelecteeLane = Sum2_WaitingTime_NonSelecteeLane + bsum*bsum;
        queue_eobm_NonSelecteeLane(1) = [];
        queue_eobm_NonSelecteeLane = [queue_eobm_NonSelecteeLane WaitingTime];
    end
end

function update_sum_variables_SelecteeLane(WaitingTime,RiskValue)
    
    global m_eobm_SelecteeLane
    global queue_eobm_SelecteeLane
    global Sumy_WaitingTime_SelecteeLane
    global Sum_WaitingTime_SelecteeLane
    global Sum2_WaitingTime_SelecteeLane
    global Sumd_WaitingTime_SelecteeLane
    global Num_CountedVisitors_SelecteeLane
    global RiskSum_SelecteeLane
    %disp(WaitingTime)
    global Test_SelecteeWaitingTime
    Test_SelecteeWaitingTime = [Test_SelecteeWaitingTime WaitingTime];
    %if length(Test_SelecteeWaitingTime) > 97
        %disp('!')
    %end
    RiskSum_SelecteeLane = RiskSum_SelecteeLane + RiskValue;
    if m_eobm_SelecteeLane < 0
        return
    end
    
    if Num_CountedVisitors_SelecteeLane <= m_eobm_SelecteeLane
        %process the first m observations---the first complete batch     
        Sumy_WaitingTime_SelecteeLane = Sumy_WaitingTime_SelecteeLane + WaitingTime;
        queue_eobm_SelecteeLane = [queue_eobm_SelecteeLane WaitingTime];
        if Num_CountedVisitors_SelecteeLane == m_eobm_SelecteeLane
            Sum_WaitingTime_SelecteeLane = Sumy_WaitingTime_SelecteeLane;
            Sum2_WaitingTime_SelecteeLane = Sum_WaitingTime_SelecteeLane*Sum_WaitingTime_SelecteeLane;
            Sumd_WaitingTime_SelecteeLane = 0;
        end
    else
        %process observations m+1 through n
        Sumd_WaitingTime_SelecteeLane = Sumd_WaitingTime_SelecteeLane + queue_eobm_SelecteeLane(1);
        Sumy_WaitingTime_SelecteeLane = Sumy_WaitingTime_SelecteeLane + WaitingTime;
        bsum = Sumy_WaitingTime_SelecteeLane - Sumd_WaitingTime_SelecteeLane;
        Sum_WaitingTime_SelecteeLane = Sum_WaitingTime_SelecteeLane + bsum;
        Sum2_WaitingTime_SelecteeLane = Sum2_WaitingTime_SelecteeLane + bsum*bsum;
        queue_eobm_SelecteeLane(1) = [];
        queue_eobm_SelecteeLane = [queue_eobm_SelecteeLane WaitingTime];
    end
end

function [ AvgWait,...
      VarAvgWait,...
      AvgWait_NonSelecteeLane,...
      VarAvgWait_NonSelecteeLane,...
      AvgWait_SelecteeLane,...
      VarAvgWait_SelecteeLane] = obm
%  purpose: to compute the obm estimate to estimate var(ybar)
%  reference:  chapter ii.4 in wheyming tina song's dissertation
%  parameter definitions
%  input
%    n:     number of observations
%    m:     batch size, an algorithm parameter (0 .lt. m .lt. n)
%    sumy:  sum of data points
%    sum:   sum of the batch sums
%    sum2:  sum of the squared batch sums
%  output
%    VarAvgWait: estimated variance of the sample mean (Average Waiting
%    Time)
%   
%  variable definitions
%    ybar:  sample average
%    sumbm: sum of the batch means
%    sum2bm:sum of the squared batch means
  
    global Num_CountedVisitors
    global Num_CountedVisitors_NonSelecteeLane
    global Num_CountedVisitors_SelecteeLane
  
    global m_eobm
    global m_eobm_NonSelecteeLane
    global m_eobm_SelecteeLane

    global Sumy_WaitingTime
    global Sum_WaitingTime
    global Sum2_WaitingTime
    global Sumy_WaitingTime_NonSelecteeLane
    global Sum_WaitingTime_NonSelecteeLane
    global Sum2_WaitingTime_NonSelecteeLane
    global Sumy_WaitingTime_SelecteeLane
    global Sum_WaitingTime_SelecteeLane
    global Sum2_WaitingTime_SelecteeLane
    
    global queue_eobm_NonSelecteeLane
    global queue_eobm_SelecteeLane
    
    AvgWait = Sumy_WaitingTime/Num_CountedVisitors;
    sumbm = Sum_WaitingTime/m_eobm;
    sum2bm = (Sum2_WaitingTime/m_eobm) / m_eobm;
    VarAvgWait = ((sum2bm - AvgWait * (2*sumbm - (Num_CountedVisitors-m_eobm+1)*AvgWait))) / ...
                (((Num_CountedVisitors-m_eobm+1.)*(Num_CountedVisitors-m_eobm))/m_eobm);
                
    AvgWait_NonSelecteeLane = 0;
    VarAvgWait_NonSelecteeLane = 0;
    if Num_CountedVisitors_NonSelecteeLane > 0
        AvgWait_NonSelecteeLane = Sumy_WaitingTime_NonSelecteeLane/Num_CountedVisitors_NonSelecteeLane;
    end
    if m_eobm_NonSelecteeLane > 0
        if m_eobm_NonSelecteeLane < Num_CountedVisitors_NonSelecteeLane
            sumbm = Sum_WaitingTime_NonSelecteeLane/m_eobm_NonSelecteeLane;
            sum2bm = (Sum2_WaitingTime_NonSelecteeLane/m_eobm_NonSelecteeLane) / m_eobm_NonSelecteeLane;
            VarAvgWait_NonSelecteeLane = ((sum2bm - AvgWait_NonSelecteeLane * (2*sumbm - (Num_CountedVisitors_NonSelecteeLane-m_eobm_NonSelecteeLane+1)*AvgWait_NonSelecteeLane))) / ...
                        (((Num_CountedVisitors_NonSelecteeLane-m_eobm_NonSelecteeLane+1.)*(Num_CountedVisitors_NonSelecteeLane-m_eobm_NonSelecteeLane))/m_eobm_NonSelecteeLane);
        elseif Num_CountedVisitors_NonSelecteeLane >= 2
            m_eobm_NonSelecteeLane = floor(Num_CountedVisitors_NonSelecteeLane/2);
            VarAvgWait_NonSelecteeLane = eobm(  Num_CountedVisitors_NonSelecteeLane,...
                                                queue_eobm_NonSelecteeLane,...
                                                m_eobm_NonSelecteeLane);
        end
    end
    AvgWait_SelecteeLane = 0;
    VarAvgWait_SelecteeLane = 0;
    if Num_CountedVisitors_SelecteeLane > 0
        AvgWait_SelecteeLane = Sumy_WaitingTime_SelecteeLane/Num_CountedVisitors_SelecteeLane;
    end    
    if m_eobm_SelecteeLane > 0
        if m_eobm_SelecteeLane < Num_CountedVisitors_SelecteeLane
            AvgWait_SelecteeLane = Sumy_WaitingTime_SelecteeLane/Num_CountedVisitors_SelecteeLane;
            sumbm = Sum_WaitingTime_SelecteeLane/m_eobm_SelecteeLane;
            sum2bm = (Sum2_WaitingTime_SelecteeLane/m_eobm_SelecteeLane) / m_eobm_SelecteeLane;
            VarAvgWait_SelecteeLane = ((sum2bm - AvgWait_SelecteeLane * (2*sumbm - (Num_CountedVisitors_SelecteeLane-m_eobm_SelecteeLane+1)*AvgWait_SelecteeLane))) / ...
                        (((Num_CountedVisitors_SelecteeLane-m_eobm_SelecteeLane+1.)*(Num_CountedVisitors_SelecteeLane-m_eobm_SelecteeLane))/m_eobm_SelecteeLane);
        elseif Num_CountedVisitors_SelecteeLane >= 2
            m_eobm_SelecteeLane = floor(Num_CountedVisitors_SelecteeLane/2);
            VarAvgWait_SelecteeLane = eobm(  Num_CountedVisitors_SelecteeLane,...
                                                queue_eobm_SelecteeLane,...
                                                m_eobm_SelecteeLane);
        end
    end
    %disp('Test')
end
