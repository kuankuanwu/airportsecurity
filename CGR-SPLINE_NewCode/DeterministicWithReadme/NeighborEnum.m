function [n_prime2,x,g] = NeighborEnum(x0)
%  purpose:         A subroutine in SPLINE method (NE). To traverse the
%                   neighbor point (S1,S2) and find the local minimum of
%                   tau.
%  input:
%   x0:             The 2*1 vector contains (S1,S2).
%  output:
%   n_prime2:       The traversed point number in this procedure. Since
%                   this procedure traverse ([S1-1,S1+1],[S2-1,S2+1])
%                   therefore if it is an interior point (i.e. the neighbor
%                   points are in feasible region) the n_prime2 is equal 9
%                   after running this procedure.
    n_prime2 = 0;
    g = Inf;
    X_best = [];
    for i = -1 : 1
		if i == 0
			for j = -1 : 1
				n_prime2 = n_prime2 + 1;
				S1 = x0(1) + i;
				S2 = x0(2) + j;
				g_cur = GetOptimTau(S1,S2,1);
				if g_cur < g
				   X_best = [S1; S2];
				   g = g_cur;
				end
			end
		else
			j = 0;
			S1 = x0(1) + i;
			S2 = x0(2) + j;
			g_cur = GetOptimTau(S1,S2,1);
			if g_cur < g
			   X_best = [S1; S2];
			   g = g_cur;
			end
		end
    end
    x = X_best;
    return
end