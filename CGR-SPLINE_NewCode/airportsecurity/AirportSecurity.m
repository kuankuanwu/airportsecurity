function [flag2, xstruct, iseed] = AirportSecurity(param, x, m, iseed)
	
% INPUT
%       param
%         param(1) = binary (yes=1/no=1) = 0  
%         param(2) = problem dimension = 4
%         param(3) = nseeds = 3
%         param(4) = nSecMeas = 1 
%         param(5) = warm-up time
%         param(6) = simulation end time (not in use)
%         param(7) = total service rate (not in use)
%         param(8) = total buffer space available (not in use)
%       x 
%         x(0) = dimension d=2
%         integer solution, ( s1, s2)
%       m = sample size
%       iseed  
%         iseed(0) = number of random-number seeds, nseeds = 3
%         iseed(1), iseed(2), iseed(3) = one seed per server
% OUTPUT
%       flag1  = 0 implies that the model parameters are feasible
%              = 1 implies that the model parameters are infeasible
%       flag2  = 0 implies that x is feasible
%              = 1 implies that x is infeasible
%       fn     = ybar (defined only if flag1 = 0 and flag2 = 0)	
%       constraint = matrix (size = 1 X nsecMeas) of estimates of 
%                   constraint functions 
%       <...>  = see oracle.m for details
%(Constraint)
global c1
global c2
global Beta1
global Beta2
global Budget
global Epsilon


xstruct=struct('x', {}, 'fn', {}, 'FnVar', {}, 'FnGrad', {}, 'FnGradCov', {}, 'constraint', {}, 'ConstraintCov', {}, 'ConstraintGrad', {}, 'ConstraintGradCov', {});
xstruct(1).x=x;

global mobs
global SecurityLevel_Global
global Constraints_Global

fn=2; %Infeasible


Upper = 0;
Lower = 0;
Upper = floor((Budget - Beta1 - c1*x(1))/c2);
Lower = floor(max((Budget - Beta2 - c1*x(1))/c2,0));

if x(2) > Upper || x(2) < Lower
  flag2 = 1;
    return  
end

Delta0 = param(6);
ErrorTol = param(7)/(2^(1/2)); %% 12/3


disp(x)
disp(iseed)
disp(m)
disp(['param(7)=' num2str(param(7)) ])
disp(['ErrorTol=' num2str(ErrorTol)])


[tau,AvgWaitingTime,Se_AvgWaitingTime,MinTau_Budget,flag2,mobs,iseed,...
    SecurityLevel_Global,Constraints_Global] = GetOptimTau_Simulation_Bisection(x(1),x(2),iseed,...
    param(5),m,...
    Delta0,ErrorTol);

if flag2 == 0
    fn = tau;
    FnVar=0;
    %P = CalcP(tau,theta);
    %constraint = [AvgWaitingTime*60 - Epsilon c1*x(1) + c2*x(2) + Beta1*P + Beta2*(1-P) - Budget];
    %ConstraintCov = [(Se_AvgWaitingTime*60)^2 0];
    constraint = [AvgWaitingTime - Epsilon];
    ConstraintCov = [(Se_AvgWaitingTime)^2 ];
    xstruct(1).fn=fn;
    xstruct(1).FnVar=FnVar;
    xstruct(1).constraint = constraint;
    xstruct(1).ConstraintCov=ConstraintCov;


end

end