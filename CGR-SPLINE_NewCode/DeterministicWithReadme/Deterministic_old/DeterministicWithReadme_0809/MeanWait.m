function [WaitingTime, PE1, PE2, PE12] = MeanWait(tau,S1,S2)
%  purpose:     To evaluate the mean waiting time given tau, S1 and S2.
%               Two methods including chain and matrix inverse are
%               implemented. To avoid the error propagation in chain
%               method, the default is set as matrix inverse.
%  input:
%   tau:            The min possible tau which satisfy the budget and mean
%                   waiting time constraint.
%   S1:             The number of servers in non-selectee lane.
%   S2:             The number of servers in selectee lane.
%  output:
%   WaitingTime:    The mean waiting time given tau, S1 and S2.
%   PE1:            The probability of the state which the visitors can only
%                   enter the non-selectee lane
%   PE2:            The probability of the state which the visitors can only
%                   enter the selectee lane
%   PE12:           The probability of the state which the visitors can
%                   enter both non-selectee lane and selectee lane.
    global MeanWait_Flag
    if ~isempty(MeanWait_Flag)
        if MeanWait_Flag == 2
            [WaitingTime, PE1, PE2, PE12] = MeanWait_Chain(tau,S1,S2);
        else
            [WaitingTime, PE1, PE2, PE12] = MeanWait_MatrixInverse(tau,S1,S2);
        end
    else
        [WaitingTime, PE1, PE2, PE12] = MeanWait_MatrixInverse(tau,S1,S2);
    end
end