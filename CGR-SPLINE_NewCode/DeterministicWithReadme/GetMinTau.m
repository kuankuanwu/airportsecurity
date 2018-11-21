function mintau = GetMinTau(S1,S2)
%  purpose: To evaluate the min possible tau which satisfy the budget
%  constraint
%  input
%   S1:         The number of servers in non-selectee lane.
%   S2:         The number of servers in selectee lane.
%  output
%   mintau:     The min possible tau which satisfy the budget. If all  
%               tau in [0,1] satisfies the budget, it returns a value less than or
%               equal to 0. If all tau in [0,1] does not satisfy the budget,
%               it returns a imaginary value.
%  variables:
%   c1:         The cost per server in non-selectee lane
%   c2:         The cost per server in selectee lane
%   Beta1:      The depreciation of the machine in non-selectee lane
%   Beta2:      The depreciation of the machine in selectee lane
%   theta:      Parameters of truncaed exponential distribution which is
%               assumed as the distriubtion of visitors' risk value
%   Budget:     Budget constraint

    global c1
    global c2
    global Beta1
    global Beta2
    global theta
    global Budget
    mintau = max(-theta*log( 1 - (1-exp(-1/theta))/(Beta1-Beta2)*(Budget - Beta2 - c1*S1 - c2*S2)),0);
    return
end