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
    global d1
    global d2
    
    global curr_s1s2
    global curr_tau
    global curr_Constraint_Time
    global curr_Constraint_Budget
    global best_s1s2
    global best_tau
    global best_Constraint_Time
    global best_Constraint_Budget
    global Awork
    global Apoints
    
    global ncalls
    global npoints
    global seed
    global Lambda 
    global mu1
    global mu2 
    
    global Tol
    global MaxIter
    global Tol_wait_fmin
    global Iter_fmin
    param=[c1,c2,Beta1,Beta2,Budget,Epsilon,theta,Lambda,mu1,mu2];
%     param=["c1",'c2','beta1','beta2','Budget','Epsilon','theta','Lambda','mu1','mu2';param];
    curr_s1s2 = [];
    curr_tau = [];
    curr_Constraint_Time = [];
    curr_Constraint_Budget = [];
    best_s1s2 = [];
    best_tau = [];
    best_Constraint_Time = [];
    best_Constraint_Budget = [];
    Awork = [];
    Apoints = [];
    
    ncalls = 0;
    npoints = 0;
    
    MaxS2_Loose = floor((Budget - Beta1)/c2);
    MaxS1_Loose = floor((Budget - Beta1)/c1);

    %RefMatrix: The record of computed optimized tau for given (S1,S2)
    %The purpose is not to recompute tau for (S1,S2) that has been
    %searched.
    global RefMatrix
    RefMatrix = -ones(MaxS1_Loose,MaxS2_Loose); 

    global OutputMatrix
    OutputMatrix = zeros(1,10);
    
    %Init = zeros(Iter,2);
    s_best = [];
    tau_best = Inf;
    OUTPUT=[];
    i = 1;
    rng('default'); %Fixed the initial random number generator
    rng(seed);
    while 1 == 1
        Tol = 1E-6/1.01^(i);    
        MaxIter=ceil(1000*1.01^(i));
        Iter_fmin=ceil(100*1.01^(i));
        Tol_wait_fmin=1E-2/1.01^(i);
        B_r=ceil(b*1.01^(i));
        %Generate an observation of S2
        while 1 == 1
            CurS1 = randi([0 MaxS1_Loose]);
            %Compute the upper and lower bounds of S1 for this S2 value
            Upper = floor((Budget - Beta1 - c1*CurS1)/c2);
            Lower = floor((Budget - Beta2 - c1*CurS1)/c2);
            Lower = max(Lower,0);
            %Generate an obervation of S1 within upper and lower bounds
            CurS2 = randi([Lower Upper]);
            if CurS2*CurS1 > 0
                break
            end
        end
        
        OutputMatrix(1,1)=CurS1;
        OutputMatrix(1,2)=CurS2;
        
        %Init(i,:) = [CurS1 CurS2];
        disp(['Iter: ' num2str(i-1) ' start: (' num2str(CurS1) ',' num2str(CurS2) ')'])
        [s, tau, n] = SPLINE([CurS1; CurS2],B_r);
        if isinf(tau)
            continue
        end
        
        OutputMatrix(1,3)=s(1);   
        OutputMatrix(1,4)=s(2);
        OutputMatrix(1,5)=tau;
        
        max_securitylevel = d1*CalcR1(tau,theta) + d2*CalcR2(tau,theta);
        OutputMatrix(1,10)=max_securitylevel;
        
        curr_s1s2 = [curr_s1s2 s];
        curr_tau = [curr_tau tau];
        Time = MeanWait(tau,s(1),s(2));
        curr_Constraint_Time = [curr_Constraint_Time (Time - Epsilon)];
        OutputMatrix(1,6)=Time;
        
        Cost = c1*s(1) + c2*s(2) + Beta1*CalcP(tau,theta) + Beta2*(1-CalcP(tau,theta));
        curr_Constraint_Budget = [curr_Constraint_Budget Cost-Budget];
        OutputMatrix(1,7)=(Cost-Budget);
       
        Awork = [Awork ncalls];
        Apoints = [Apoints npoints];
        OutputMatrix(1,9)=ncalls;
        OutputMatrix(1,8)=npoints;
        
        
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
        OUTPUT=[OUTPUT;OutputMatrix];
        i = i + 1;
%         xlswrite(strcat('Deter_c-',num2str(c2/c1),'_beta-',num2str(Beta2/Beta1),'_mu-',num2str(mu1/mu2),'beta-',num2str(Beta2/Beta1),'_Lambda-',num2str(Lambda),'.xls'),OutputMatrix,'¤u§@ªí1',strcat('A',int2str(i)));
        if i > Iter
            break;
        end
    end
    
%     xlswrite(strcat('Deter_c1-',num2str(c1),'c2-',num2str(c2),'_beta1-',num2str(Beta1),'beta2-',num2str(Beta2),'.xls'),OutputMatrix,'result');
%     xlswrite(strcat('Deter_c1-',num2str(c1),'c2-',num2str(c2),'_beta1-',num2str(Beta1),'beta2-',num2str(Beta2),'.xls'),param,'param');
  
    x_best =  OUTPUT;
 

    
    
end