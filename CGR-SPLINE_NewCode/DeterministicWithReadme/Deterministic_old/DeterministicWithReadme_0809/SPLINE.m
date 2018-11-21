function [x, g, n] = SPLINE(x0,bk)
%  purpose: To search the min possible tau which satisfies both the mean waiting time
%           constraint and budget constraint by estimating the gradient and
%           enumerating the neighbor points.
%           if the results from gradient method is the same as the results
%           from neighbor points enumeration, this method will return this
%           results and terminate.
%  input:
%   x0:     The 2*1 vector contains (S1,S2).
%   bk:     Terminate condition.
%           The limit of traversed points (S1,S2). If the optimization
%           procedure evalute the min possible tau more than bk times, this
%           procedure will terminate and return the current solution.
%  output:
%   x:      The (S1,S2) point which can produce min possible tau 
%           satisfying both the mean waiting time constraint and budget constraint
%   g:      The min possible tau satisfies both constraints given (S1,S2).
    global RefMatrix
    
    n = 0;
    if isempty(GetOptimTau(x0(1),x0(2),1))
        disp('Not Feasible x0: Try other feasible initial point')
        x = x0;
        g = Inf;
        return
    end
    
    x_new = x0;
    while 1 == 1
        %g_x_new = GetOptimTau2(x0(1),x0(2),1);
        [n_prime,x_old,g_x_old] = SPLI(x_new,bk);
        disp(['xold: ' num2str(x_old(1)) ',' num2str(x_old(2)) ])
        [n_prime2,x_new,g_x_new] = NeighborEnum(x_old);
        disp(['xnew: ' num2str(x_new(1)) ','  num2str(x_new(2)) ])
        n = n + n_prime + n_prime2;
        if g_x_new == g_x_old
%           if the results from gradient method is the same as the results
%           from neighbor points enumeration, this method will return this
%           results and terminate.
            x = x_new;
            g = GetOptimTau(x(1),x(2),1);
            return
        end
        if n > bk
            x = x_new;
            g = GetOptimTau(x(1),x(2),1);
            return
        end
    end
    
    g = GetOptimTau(x(1),x(2),1);
end