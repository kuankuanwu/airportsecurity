function [ub_guess, AvgWaitingTime_Guess, Se_AvgWaitingTime_Guess] = UpperBoundGuess_Simulation(S1,S2,Epsilon,iseed0,Num_Warmup,Num_Simulation)
%  purpose: To compute the upper bound (i.e. the tau with the corresponded Avg waiting time < epsilon) 
%           for regular falsi method which finds the minimum tau satisfies the constraint of Avg waiting time.
%           This method assumes that the tau which makes the ratio of assigned visitors in two lanes equals to
%           the ratio of service capabilities in two lanes can provide the
%           highest efficiency to serve the visitors and hence the putative lowest Avg waiting time can be obtained.
%
%           The putative lowest Avg waiting time usually close to the actual lowest
%           Avg waiting time and satisfy the constraint of Avg waiting time. If
%           this putative value cannot satisfy the constraint of Avg waiting time,
%           the fminsearch() should be applied to find the actual Avg waiting time.
%
%  input
%   S1:         The number of servers in non-selectee lane
%   S2:         The number of servers in selectee lane
%   Epsilon:    The constraint of Avg waiting time.
%  output
%    ub:        The upper bound for regular falsi method which finds the minimum
%               tau satisfies the constraint of Avg waiting time.
%    flag:      The indicator of the validity of ub. Valid if flag > 0
%               otherwise the ub is invalid which Avgs the putative Avg waiting 
%               time cannot satisfy the constraint.
%  variables:
%   mu1:        The service rate per server in non-selectee lane
%               (visitors/hr)
%   mu2:        The service rate per server in selectee lane
%               (visitors/hr)
%   theta:      Parameters of truncaed exponential distribution which is
%               assumed as the distriubtion of visitors' risk value


    global mu1
    global mu2
    global theta
    p_guess = mu1*S1/(mu1*S1 + mu2*S2);
    tau_guess = -theta*log(1 - p_guess*(1 - exp(-1/theta)));
    ub_guess = tau_guess;
    iseed_vector = iseed0(1:4);
    [AvgWaitingTime_Guess,Se_AvgWaitingTime_Guess] = Simulation_AirportModel(tau_guess,S1,S2,iseed_vector,Num_Warmup,Num_Simulation);
end