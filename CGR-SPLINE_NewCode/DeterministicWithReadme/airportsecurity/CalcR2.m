function R2 = CalcR2(tau, theta)
    R2 = (exp(-1/theta) - tau*exp(-tau/theta) + theta*exp(-1/theta) - theta*exp(-tau/theta));
    R2 = R2/(exp(-1/theta) + theta*exp(-1/theta) - theta);
end