function [tau_minWT,lb,ub, tau_minWT_func, lb_func, ub_func, iseed, ncalls_search] = FindMinAvgWait(Epsilon,S1,S2,tau_budget,AvgWaitingTime_tau_budget, tau_guess, AvgWaitingTime_tau_guess,iseedk,...
                Num_Warmup,Num_Simulation,...
                Tol,MaxIter)
% Purpose: Finding bounding intervals of the smallest root            
            
            
global log_test
ncalls_search = 0;
[AvgWaitingTime_0,Se_AvgWaitingTime_0,iseed] = Simulation_AirportModel(0,S1,S2,iseedk,Num_Warmup,Num_Simulation);
[AvgWaitingTime_1,Se_AvgWaitingTime_1,iseed] = Simulation_AirportModel(1,S1,S2,iseedk,Num_Warmup,Num_Simulation);
ncalls_search = ncalls_search + Num_Simulation*2 + Num_Warmup*2;
x = [0; tau_budget; tau_guess; 1];
y = [AvgWaitingTime_0; AvgWaitingTime_tau_budget; AvgWaitingTime_tau_guess; AvgWaitingTime_1];
CurIter = 0;

lb = tau_budget;
lb_func = AvgWaitingTime_tau_budget;
old_x = Inf;


while 1 == 1
    disp('====Start====')
    fprintf(log_test,'====Start====\n');
    %Select three point from four initial point set:
    %Sort by tau
    Tmp = sortrows([x,y]);
    x = Tmp(:,1);
    y = Tmp(:,2);
    disp([x y])
    fprintf(log_test,'(x1,y1) = (%f,%f)\n',x(1),y(1));
    fprintf(log_test,'(x2,y2) = (%f,%f)\n',x(2),y(2));
    fprintf(log_test,'(x3,y3) = (%f,%f)\n',x(3),y(3));
    fprintf(log_test,'(x4,y4) = (%f,%f)\n',x(4),y(4));
    %1. Select the min
    MinIdx = 1;
    Flag_OnlyThreePoints = 0;
    for i = 2 : length(y)
        if y(i) < y(MinIdx)
            MinIdx = i;
        end
        if x(i) - x(i-1) < Tol
            Flag_OnlyThreePoints = 1;
            x(i) = [];
            y(i) = [];
            break;
        end
    end
    if Flag_OnlyThreePoints == 0
        if MinIdx <= 2
            x = x(1:3);
            y = y(1:3);
        else
            x = x(2:4);
            y = y(2:4);
        end
    end
    disp('====Get Three Points====')
    fprintf(log_test,'====Get Three Points====\n');
    fprintf(log_test,'(x1,y1) = (%f,%f)\n',x(1),y(1));
    fprintf(log_test,'(x2,y2) = (%f,%f)\n',x(2),y(2));
    fprintf(log_test,'(x3,y3) = (%f,%f)\n',x(3),y(3));
    disp([x y])
    % Find the quadratic function ax^2+bx+c, going through these three points
%     A = zeros(3,3);
%     for i = 1 : 3
%         A(i,:) = [x(i)^2 x(i) 1];
%     end
%     coef = linsolve(A,y);
%     new_x = -coef(2)/(2*coef(1));
    a = (((y(3)-y(2))/(x(3)-x(2))) - ((y(3)-y(1))/(x(3)-x(1)))) / (x(2) - x(1));
    b = ((y(3)-y(2))/(x(3)-x(2))) - (x(3)+x(2))*a;
    new_x = -b/(2*a);
    disp('====Get New Points====')
    if new_x > 1 || new_x < 0
        disp('ERROR')
        fprintf(log_test,'Error: (x_new) = (%f)\n',new_x);
    end
    [new_y,new_y_se,iseed] = Simulation_AirportModel(new_x,S1,S2,iseedk,Num_Warmup,Num_Simulation);
    ncalls_search = ncalls_search + Num_Simulation + Num_Warmup;
    disp([new_x new_y])
    fprintf(log_test,'New Points: (x_new,y_new) = (%f,%f)\n',new_x,new_y);
    if new_y < Epsilon && new_x > tau_budget
        ub = new_x;
        ub_func = new_y;
        tau_minWT = new_x;
        tau_minWT_func = new_y;
        for i = 1 : 3
            if x(i) > lb && x(i) < ub
                lb = x(i);
                lb_func = y(i);
            end
        end
        disp('Success')
        fprintf(log_test,'Done!\n');
        return
    else
        
        if (max([x; new_x]) - min([x; new_x])) < Tol || abs(old_x - new_x) < Tol
            %Converged
            Tmp_x = [x; new_x];
            Tmp_y = [y; new_y];
            Idx = find(Tmp_y == min(Tmp_y));
            Idx = Idx(1);
            ub = [];
            ub_func = [];
            tau_minWT = Tmp_x(Idx);
            tau_minWT_func = Tmp_y(Idx);
            fprintf(log_test,'Return Min Point: (x_min, y_min) = (%f,%f)\n', tau_minWT, tau_minWT_func);
            return
        end
        
        
    end
    
    x = [x; new_x];
    y = [y; new_y];
    old_x = new_x;
    disp(['EndIter ' num2str(CurIter)])
    fprintf(log_test,'End iter %d\n', CurIter);
    CurIter = CurIter + 1;

    if CurIter > MaxIter
        %Cannot get lower avgWT
        Idx = find(y == min(y));
        Idx = Idx(1);
        ub = [];
        ub_func = [];
        tau_minWT = x(Idx);
        tau_minWT_func = y(Idx);
        fprintf(log_test,'Reach Maximum Iteration Num: (x_min, y_min) = (%f,%f)\n', tau_minWT, tau_minWT_func);
        return
    end
end
end
    
    