    %  Driver program of airport security system optimization
%  purpose:     To optimize the airport security system given the
%               corresponded parameters including the cost, constraint and budget.
%               After the parameters were set, the GlobalOptimization() is
%               applied to find the best parameters (S1,S2,tau) using the SPLINE method
%               with different initial points.
%  input file
%   Input.txt:  Contains the parameters of the systems
%  output file:
%   Output.txt: Log file during optimization. Finally it will record the
%               minimum tau and the maximum security level.



global theta

global NInitial
global b
global d1
global d2

global fout

global curr_s1s2
global curr_tau
global curr_Constraint_Time
global curr_Constraint_Budget
global best_s1s2
global best_tau
global best_Constraint_Time
global best_Constraint_Budget
global chainflag
global N
for i=1:1
tic
InputFileName = ['input',num2str(i),'.txt'];
InitializeAndReadInput(InputFileName);
chainflag =0;
%----------------Start------------------------

if chainflag==1
    N=200;
end
%b = 10000; %Maximum number of spline iterates, used in SPLINE.m

x_best = GlobalOptimization(NInitial,b);
best_S1 = x_best(1);    
best_S2 = x_best(2);
best_tau = x_best(3);
time=toc;
max_securitylevel = d1*CalcR1(best_tau,theta) + d2*CalcR2(best_tau,theta);
save(strcat(InputFileName, 'Case33_vars.mat'), 'x_best', 'best_S1', 'best_S2', 'best_tau','time');
disp('')
disp('=========================================================================')
disp(['Global optimal point (S1,S2,tau) = (' num2str(best_S1) ',' num2str(best_S2) ',' num2str(best_tau) ')'])
disp(['The optimal security level = ' num2str(max_securitylevel)])

fprintf(fout,'==========================================================')
fprintf(fout,'Global optimal point (S1,S2,tau) = (%d,%d,%f)\n',best_S1,best_S2,best_tau)
fprintf(fout,'The optimal security level = %e\n',max_securitylevel)
fclose(fout)

end

