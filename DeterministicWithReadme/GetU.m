function U = GetU(N,S1)
    U = zeros(S1+N+1,S1+N+1);
    U(1:S1+N,1:S1+N) = eye(S1+N);
    return
end