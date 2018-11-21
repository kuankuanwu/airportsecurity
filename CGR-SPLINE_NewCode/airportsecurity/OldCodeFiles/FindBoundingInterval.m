function [lb,ub,flag, iseed, ncalls_search] = FindBoundingInterval(Epsilon,S1,S2,tau_budget,meanWaitingTime_tau_budget, tau_guess, meanWaitingTime_tau_guess,iseedk,...
                Num_Warmup,Num_Simulation,...
                Delta0,Tol,MaxIter)
% Purpose: Finding bounding intervals of the smallest root            
            
            
global log_test
ncalls_search = 0;
[meanWaitingTime_0,Se_AvgWaitingTime_0,iseed] = Simulation_AirportModel(0,S1,S2,iseedk,Num_Warmup,Num_Simulation);
[meanWaitingTime_1,Se_AvgWaitingTime_1,iseed] = Simulation_AirportModel(1,S1,S2,iseedk,Num_Warmup,Num_Simulation);
ncalls_search = ncalls_search + Num_Simulation*2;
x = [0; tau_budget; tau_guess; 1];
y = [meanWaitingTime_0; meanWaitingTime_tau_budget; meanWaitingTime_tau_guess; meanWaitingTime_1];
CurIter = 0;

while 1 == 1
    disp('====Start====')
    disp([x y])
    %Select three point from four initial point set:
    %Sort by tau
    Tmp = sortrows([x,y]);
    x = Tmp(:,1);
    y = Tmp(:,2);
    %1. Select the min
    MinIdx = 1;
    for i = 2 : length(y)
        if y(i) < y(MinIdx)
            MinIdx = i;
        end
    end
    if MinIdx <= 2
        x = x(1:3);
        y = y(1:3);
    else
        x = x(2:4);
        y = y(2:4);
    end
    disp('====Get Three Points====')
    disp([x y])
    % Find the quadratic function ax^2+bx+c, going through these three points
    A = zeros(3,3);
    for i = 1 : 3
        A(i,:) = [x(i)^2 x(i) 1];
    end
    coef = linsolve(A,y);
    new_x = -coef(2)/(2*coef(1));
    disp('====Get New Points====')
    if new_x > 1 || new_x < 0
        disp('ERROR')
    end
    [new_y,new_y_se,iseed] = Simulation_AirportModel(new_x,S1,S2,iseedk,Num_Warmup,Num_Simulation);
    ncalls_search = ncalls_search + Num_Simulation;
    disp([new_x new_y])
    if new_y*60 < Epsilon && new_x > tau_budget
        ub = new_x;
        lb = tau_budget;
        for i = 1 : 3
            if x(i) > tau_budget && x(i) < ub
                lb = x(i);
            end
        end
        flag = 1;
        disp('Success')
        return
    else
        if new_y > min(y)
            %Cannot get lower avgWT
            %lb = [];
            %ub = [];
            %flag = -1;
            %return
        else
            if min(y) - new_y < Tol
                %Cannot get lower avgWT
                lb = [];
                ub = [];
                flag = -1;
                return
            end
        end

    end
    
    x = [x; new_x];
    y = [y; new_y];
    disp(['EndIter ' num2str(CurIter)])
    CurIter = CurIter + 1;

    if CurIter > MaxIter
        %Cannot get lower avgWT
        lb = [];
        ub = [];
        flag = -1;
        return
    end
end
end
    
    