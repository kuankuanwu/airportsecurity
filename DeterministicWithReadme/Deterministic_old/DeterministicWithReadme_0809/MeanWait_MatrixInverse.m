function [WaitingTime, PE1, PE2, PE12] = MeanWait_MatrixInverse(tau,S1,S2)
%  purpose:         To evaluate the mean waiting time given tau, S1 and S2
%                   using matrix inverse method.
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
    global OldQN2
    global N
    global Lambda
    global mu1
    global mu2
    global theta
    
    global ncalls
    
    PE1 = -1;
    PE2 = -1;
    PE12 = -1;
    if tau < 0 
        WaitingTime = MeanWait_MatrixInverse(0,S1,S2) + abs(tau);
        return
    end
    if tau > 1
        WaitingTime = MeanWait_MatrixInverse(1,S1,S2) + abs(tau);
        return
    end
    if S1 == 0 && tau ~= 0 %Send visitor to non-selected lane but no server
        WaitingTime = Inf;
        return
    end
    if S2 == 0 && tau ~= 1 %Send visitor to selected lane but no server
        WaitingTime = Inf;
        return
    end
    
    ncalls = ncalls + 1;
    
    p = (-exp(-tau/theta) + 1)/(1 - exp(-1/theta));
    if S1 == 0 && tau == 0
        QN2 = zeros(1*(S2+N+1),1*(S2+N+1));
        A{S2+N+1,S2+N+1} = [];
        Lambda2 = Lambda*(1-p);
        for i = 0 : (S2+N)
            base = i*1;
            if i < (S2 + N)
                ArrivalRate = Lambda2;
                base_right = base + 1;
                QN2(base + (1:1), base_right + (1:1)) = ArrivalRate;
                A{i +1,i+1 +1} = QN2(base + (1:1), base_right + (1:1));
                QN2(base + (1:1), base + (1:1)) = -ArrivalRate;
                A{i +1,i +1} = QN2(base + (1:1), base + (1:1));
            else
            end
            if i > 0
                if i < S2
                    ProcessRate = mu2 * i;
                else
                    ProcessRate = mu2 * S2;
                end
                base_left = base - 1;
                QN2(base + (1:1), base_left + (1:1)) = ProcessRate;
                A{i +1,i-1 +1} = QN2(base + (1:1), base_left + (1:1));
                
                QN2(base + (1:1), base + (1:1)) = QN2(base + (1:1), base + (1:1)) - ProcessRate;
                A{i +1,i +1} = QN2(base + (1:1), base + (1:1));
            end
        end
    elseif S2 == 0 && tau == 1
        QN2 = zeros((S1+N+1)*1,(S1+N+1)*1);
        A{S1+N+1,S1+N+1} = [];
        Lambda1 = Lambda*p;
        for i = 0 : (S1+N)
            base = i*1;
            if i < (S1 + N)
                ArrivalRate = Lambda1;
                base_right = base + 1;
                QN2(base + (1:1), base_right + (1:1)) = ArrivalRate;
                A{i +1,i+1 +1} = QN2(base + (1:1), base_right + (1:1));
                QN2(base + (1:1), base + (1:1)) = -ArrivalRate;
                A{i +1,i +1} = QN2(base + (1:1), base + (1:1));
            else
            end
            if i > 0
                if i < S1
                    ProcessRate = mu1 * i;
                else
                    ProcessRate = mu1 * S2;
                end
                base_left = base - 1;
                QN2(base + (1:1), base_left + (1:1)) = ProcessRate;
                A{i +1,i-1 +1} = QN2(base + (1:1), base_left + (1:1));
                
                QN2(base + (1:1), base + (1:1)) = QN2(base + (1:1), base + (1:1)) - ProcessRate;
                A{i +1,i +1} = QN2(base + (1:1), base + (1:1));
            end
        end
    else
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
    end
    
    
    %QN2
    %Flag = sum(QN2') ~= 1;
    %QN2 = QN2(Flag,Flag);
    
    
    Tmp = QN2';
    IdxA = (sum(abs(Tmp)) ~= 0);
    IdxB = (sum(abs(Tmp')) ~= 0);
    Tmp = Tmp(IdxB, IdxA);
    Tmp = [Tmp ; ones(1,size(Tmp,2))];
    %Tmp(1,:) = [];
    TmpV = zeros(size(Tmp,1), 1);
    TmpV(end) = 1;
    Tmp = sparse(Tmp);
    TmpV = sparse(TmpV);
    
    if S1 == 0 && tau == 0
        TmpDistribution = zeros(1,1*(S2+N+1));
        TmpDistribution(IdxA) = full(Tmp\TmpV);
        TmpDistribution = reshape(TmpDistribution,1,(S2+N+1));
        Distribution = zeros(S1+N+1,S2+N+1);
        Distribution(1,:) = TmpDistribution;
    elseif S2 == 0 && tau == 1
        TmpDistribution = zeros(1,(S1+N+1)*1);
        TmpDistribution(IdxA) = full(Tmp\TmpV);
        TmpDistribution = reshape(TmpDistribution,(S1+N+1),1);
        Distribution = zeros(S1+N+1,S2+N+1);
        Distribution(:,1) = TmpDistribution;
    else
        Distribution = zeros(1,(S1+N+1)*(S2+N+1));
        Distribution(IdxA) = full(Tmp\TmpV);
        Distribution = reshape(Distribution,(S1+N+1),(S2+N+1));
    end
    
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
    
    
    PF = 0;
    for i = 0 : (size(Distribution,1)-1)
        j = N+S1+S2-i;
        if (j + 1) <= size(Distribution,2)
            PF = PF + Distribution(i+1,j+1);
        end
    end
    %PE1: P{The second station has queue length = N, but there are free servers in the first station}
    PE1 = sum(Distribution((0  +1): (S1-1  +1),end));
    %PE2: P{The first station has queue length = N, but there are free
    %servers in the second station}
    PE2 = sum(Distribution(end,(0  +1): (S2-1  +1)));
    
    PE12 = 1 - PF - PE1 - PE2 ;
    
    %Do not use tau == 0 or tau == 1
    if p == 0
        W1 = 0;
        W2 = L2/(Lambda*(1-p)*(PE2 + PE12));
    elseif p == 1
        W1 = L1/(Lambda*p*(PE1 + PE12));
        W2 = 0;
    else
        W1 = L1/(Lambda*p*(PE1 + PE12));
        W2 = L2/(Lambda*(1-p)*(PE2 + PE12));
    end
    
    Tmp1 = p*(PE1 + PE12);
    Tmp2 = (1-p)*(PE2 + PE12);
    final_p = Tmp1/(Tmp1+Tmp2);
    
    W_Test = final_p*W1 + (1-final_p)*W2;

    WaitingTime = W_Test;
    
    return