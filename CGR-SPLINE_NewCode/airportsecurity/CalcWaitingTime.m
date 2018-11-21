function WaitingTime = CalcWaitingTime(S1,S2,tau)
    global OldQN2
    %S1 = 3;
    %S2 = 4;
    %Lambda = 18;
    %mu1 = 8;
    %mu2 = 6;
    %tau = 0.05;
    %theta = 0.0625;
    N = 100;
    %S1 = 3;
    %S2 = 4;
    Lambda = 18;
    mu1 = 8;
    mu2 = 6;
    %tau = 0.05;
    theta = 0.0625;
    
    QN2 = zeros((S1+N+1)*(S2+N+1),(S1+N+1)*(S2+N+1));
    A{S2+N+1,S2+N+1} = [];
    for i = 0 : (S2+N)
        base = i*(S1+N+1);
        QN2(base + (1:(S1+N+1)), base + (1:(S1+N+1))) = GetA(i, i, N, S1, S2, Lambda, mu1, mu2, tau, theta);
        A{i +1,i +1} = QN2(base + (1:(S1+N+1)), base + (1:(S1+N+1)));
        if i > 0
            base_left = base - (S1 + N + 1);
            QN2(base + (1:(S1+N+1)), base_left + (1:(S1+N+1))) = GetA(i, i - 1, N, S1, S2, Lambda, mu1, mu2, tau, theta);
            A{i +1,i-1 +1} = QN2(base + (1:(S1+N+1)), base_left + (1:(S1+N+1)));
        end
        if i < (S2 + N)
            base_right = base + (S1 + N + 1);
            QN2(base + (1:(S1+N+1)), base_right + (1:(S1+N+1))) = GetA(i, i + 1, N, S1, S2, Lambda, mu1, mu2, tau, theta);
            A{i +1,i+1 +1} = QN2(base + (1:(S1+N+1)), base_right + (1:(S1+N+1)));
        end
    end
    
    %QN2
    %Flag = sum(QN2') ~= 1;
    %QN2 = QN2(Flag,Flag);
    Tmp = QN2';
    Tmp = [Tmp ; ones(1,size(Tmp,2))];
    Tmp(1,:) = [];
    TmpV = zeros(size(Tmp,1), 1);
    TmpV(end) = 1;
    
    Tmp_Compact = sparse(Tmp);
    TmpV_Compact = sparse(TmpV);
    Distribution = Tmp\TmpV;
    Distribution = reshape(Distribution,(S1+N+1),(S2+N+1));
    OldQN2 = QN2;
    
    
    pi_2 = sum(Distribution);
    
    W = 0;
    for j = 0 : (length(pi_2)-1) %#Visitor in Selectee
        if sum(Distribution(:,j+1)) ~= 0
            %pi_1 = Distribution(:,j+1)/sum(Distribution(:,j+1));
            pi_1 = Distribution(:,j+1);
        else
            pi_1 = Distribution(:,j+1);
        end
        for i = 0 : (length(pi_1)-1)
            W = W + i*(pi_1(i+1));
        end
    end
    for i = 0 : (length(pi_2)-1)
        W = W + i*pi_2(i+1);
    end
    W = W/Lambda;
    
    L1 = 0;
    for j = 0 : (length(pi_2)-1)
        if sum(Distribution(:,j+1)) ~= 0
            %pi_1 = Distribution(:,j+1)/sum(Distribution(:,j+1));
            pi_1 = Distribution(:,j+1);
        else
            pi_1 = Distribution(:,j+1);
        end
        for i = 0 : (length(pi_1)-1)
            L1 = L1 + i * (pi_1(i+1));
        end
    end
    
    L2 = 0;
    for i = 0 : (length(pi_2)-1)
        L2 = L2 + i*(pi_2(i+1));
    end
    
    p = (-exp(-tau/theta) + 1)/(1 - exp(-1/theta));
    W1 = L1/(Lambda*p);
    W2 = L2/(Lambda*(1-p));
    
    W_Test = p*W1 + (1-p)*W2;
    WaitingTime = W_Test;
    return