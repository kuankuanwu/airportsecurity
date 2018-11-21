function R1 = CalcR1(tau, theta)
    R1 = (tau*exp(-tau/theta) + theta*exp(-tau/theta) - theta);
    R1 = R1/(exp(-1/theta) + theta*exp(-1/theta) - theta);
end