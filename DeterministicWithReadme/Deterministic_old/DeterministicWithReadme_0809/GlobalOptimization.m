function x_best = GlobalOptimization(Iter, b)
%  purpose: To find the best parameters (S1,S2,tau) using the SPLINE method
%           with different initial points. The initial points (S1,S2) is
%           generated randomly in uniform distribution and using the SPLINE
%           to find (S1',S2') which can generate the minimal tau.
%  input
%   Iter:   Terminate condition. 
%           The limit of the number of initial points (S1,S2). Return the
%           current solution when the number of initial points is reach the
%           limit.
%   b:      Terminate condition.
%           The limit of traversed points (S1,S2). If the optimization
%           procedure evalute the min possible tau more than b times, this
%           procedure will terminate and return the current solution.
%  output
%   x_best: A 3*1 vector. The putative global optimum point (S1,S2,tau)

    global c1
    global c2
    global Beta1
    global Beta2
    global Budget
    global Epsilon
    global theta
 
    
    global curr_s1s2
    global curr_tau
    global curr_Constraint_Time
    global curr_Constraint_Budget
    global best_s1s2
    global best_tau
    global best_Constraint_Time
    global best_Constraint_Budget
    global Awork
    global evaluation
    global ncalls
    
    curr_s1s2 = [];
    curr_tau = [];
    curr_Constraint_Time = [];
    curr_Constraint_Budget = [];
    best_s1s2 = [];
    best_tau = [];
    best_Constraint_Time = [];
    best_Constraint_Budget = [];
    Awork = [];
    
    ncalls = 0;
    evaluation =0;
    
    MaxS2_Loose = ceil((Budget - Beta1)/c2);
    MaxS1_Loose = ceil((Budget - Beta1)/c1);

    %RefMatrix: The record of computed optimized tau for given (S1,S2)
    %The purpose is not to recompute tau for (S1,S2) that has been
    %searched.
    global RefMatrix
    RefMatrix = -ones(MaxS1_Loose,MaxS2_Loose); 
    
    global OutputMatrix
    OutputMatrix = zeros(10,10);
    
    %Init = zeros(Iter,2);
    s_best = [];
    tau_best = Inf;
    
    rng(3)
    
    i = 1;
    while 1 == 1
        %Generate an observation of S2
        while 1 == 1
            CurS2 = randi([0 MaxS2_Loose]);
            %Compute the upper and lower bounds of S1 for this S2 value
            Upper = floor((Budget - Beta1 - c1*CurS2)/c2);
            Lower = floor((Budget - Beta2 - c1*CurS2)/c2);
            Lower = max(Lower,0);
            %Generate an obervation of S1 within upper and lower bounds
            CurS1 = randi([Lower Upper]);
            if CurS2*CurS1 > 0
                break
            end
        end
        
        OutputMatrix(i,1)=CurS1;
        OutputMatrix(i,2)=CurS2;
        %Init(i,:) = [CurS1 CurS2];
        disp(['Iter: ' num2str(i) ' start: (' num2str(CurS1) ',' num2str(CurS2) ')'])
        [s, tau, n] = SPLINE([CurS1; CurS2],b);
        if isinf(tau)
            continue
        end
        
        OutputMatrix(i,3)=s(1);   
        OutputMatrix(i,4)=s(2);
        OutputMatrix(i,5)=tau;
        
        curr_s1s2 = [curr_s1s2 s];
        curr_tau = [curr_tau tau];
        Time = MeanWait(tau,s(1),s(2));
        curr_Constraint_Time = [curr_Constraint_Time (Time - Epsilon/60)];
        OutputMatrix(i,6)=(Time - Epsilon/60);
        
        Cost = c1*s(1) + c2*s(2) + Beta1*CalcP(tau,theta) + Beta2*(1-CalcP(tau,theta));
        curr_Constraint_Budget = [curr_Constraint_Budget Cost-Budget];
        OutputMatrix(i,7)=(Cost-Budget);
        
        Awork = [Awork ncalls];
        OutputMatrix(i,9)=ncalls;
        OutputMatrix(i,8)=evaluation;
        
        if tau < tau_best
            best_s1s2 = [best_s1s2 s];
            best_tau = [best_tau tau];
            best_Constraint_Time = [best_Constraint_Time (Time - Epsilon/60)];
            best_Constraint_Budget = [best_Constraint_Budget Cost-Budget];
            s_best = s;
            tau_best = tau;
        else
            best_s1s2 = [best_s1s2 best_s1s2(:,end)];
            best_tau = [best_tau best_tau(end)];
            best_Constraint_Time = [best_Constraint_Time best_Constraint_Time(end)];
            best_Constraint_Budget = [best_Constraint_Budget best_Constraint_Budget(end)];
        end
        disp(['Iter: ' num2str(i) ' is done!(' num2str(s(1)) ',' num2str(s(2)) ',' num2str(tau) ')' ])
        
        
        
        
        i = i + 1;
        if i > Iter
            break;
        end
    end
     
    xlswrite('de.xls',OutputMatrix);
    x_best = [s_best ; tau_best];
end