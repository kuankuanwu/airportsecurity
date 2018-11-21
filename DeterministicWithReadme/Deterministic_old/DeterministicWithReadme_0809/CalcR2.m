function R2 = CalcR2(tau, theta)
% purpose: to compute the risk ratio in selectee lane
% parameter definitions
% input
%     tau:    The threshold to decide the visitors go to non-selectee lane
%             or selectee lane. If the risk value of the visitor <= tau, this
%             visitor would be assigned to non-selectee lane otherwise to selectee lane.
%     theta:  Parameters of truncaed exponential distribution which is
%             assumed as the distriubtion of visitors' risk value
% output
%     R2:     The risk ratio in selectee lane

    R2 = (exp(-1/theta) - tau*exp(-tau/theta) + theta*exp(-1/theta) - theta*exp(-tau/theta));
    R2 = R2/(exp(-1/theta) + theta*exp(-1/theta) - theta);
end