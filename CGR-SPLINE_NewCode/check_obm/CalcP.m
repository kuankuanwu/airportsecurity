function P = CalcP(tau, theta)
    P = -exp(-tau/theta) + 1;
    P = P/(1-exp(-1/theta));
end