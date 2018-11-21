function [g_interpolate, gamma] = PLI(x)
%  purpose: A subroutine in SPLINE method (PLI). To estimate the
%           gradient and average minimal tau satisfy the mean waiting time
%           constraint and budget constraint.
%  input:
%   x:      The 2*1 vector contains perturbed (S1,S2). The perturbed S1 and S2 should
%           not be integers.
%  output:
%   g_interpolate:  Average minimal tau satisfy the mean waiting time
%                   constraint and budget constraint.
%   gamma:          The gradient ( delta(S1,S2)/delta(tau))
    x_est = zeros(length(x)+1, length(x));
    x_0 = floor(x);
    z = x - x_0;
    z = sort(z,'descend');
    
    Tmp = [(1 : length(x))' x z];
    Tmp = sortrows(Tmp,-3); %Sort the third row (z) by descend order
    w = zeros(length(x)+1,1);
    
    for i = 0 : length(x)
        if i == 0
            x_est(1,:) = x_0;
            Cur_x_est = x_0;
            w(1) = 1 - Tmp(1,3);
        else
            Cur_x_est(Tmp(i),1) = Cur_x_est(Tmp(i),1) + 1;
            x_est(i+1,:) = Cur_x_est;
            if i < length(x)
                w(i+1) = Tmp(i,3) - Tmp(i+1,3); 
            else
                w(length(x)+1) = Tmp(i,3); 
            end
        end
    end
    
    n = 0;
    t = 0;
    g_interpolate = 0;
    gamma = zeros(length(x),1);
    g_record = zeros(length(x)+1,1);
    for i = 0 : length(x)
        %Calculate min tau of this point :)
        S1 = x_est(i+1,1);
        S2 = x_est(i+1,2);
        tau = GetOptimTau(S1,S2,1);
        if isempty(tau)
            tau = Inf;
        end
        g_record(i+1) = tau;
        if tau < Inf %Feasible
            n = n+1;
            t = t + w(i+1);
            g_interpolate = g_interpolate + w(i+1)*tau;
        end
        
        if i > 0
            gamma(Tmp(i,1)) = tau - g_record(i);
        end
    end
    
    if t > 0
        g_interpolate = g_interpolate / t;
    else
        g_interpolate = Inf;
    end
    return
end