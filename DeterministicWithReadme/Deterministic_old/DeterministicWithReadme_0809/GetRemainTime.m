function RemainTime = GetRemainTime(x,S1,S2,Epsilon)
%  purpose: To evaluate the difference between mean waiting time and th
%  mean waiting time constraint. This function is only designed for
%  fminsearch() to solve the min possible tau which satisfy the mean
%  waiting time constraint
%  input
%   x:              tau value
%   S1:             The number of servers in non-selectee lane.
%   S2:             The number of servers in selectee lane.
%   Epsilon:        Mean waiting time constraint
%  output
%   RemainTime:     The difference between mean waiting time and the mean waiting time constraint
%                   (MeanWaitingTime(tau,S1,S2) - Epsilon)

    RemainTime = MeanWait_MatrixInverse(x,S1,S2)*60 - Epsilon;
    %disp(RemainTime)
    %disp(x)
end