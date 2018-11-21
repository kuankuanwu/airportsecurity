function [WaitingTime, PE1, PE2, PE12] = MeanWait_Chain(tau,S1,S2)
%  purpose:         To evaluate the mean waiting time given tau, S1 and S2
%                   using chain method.
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
    global N
    global Lambda
    global mu1
    global mu2
    global theta
    
    %tic
    if tau < 0 
        WaitingTime = MeanWait_Chain(S1,S2,0) + abs(tau);
        return
    end
    if tau > 1
        WaitingTime = MeanWait_Chain(S1,S2,1) + abs(tau);
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
    
    ReverseFlag = 0;
    if S2 == 0 && tau == 1
        %Reverse: Technical Reason
        tmp = mu2;
        mu2 = mu1;
        mu1 = tmp;
        S2 = S1;
        S1 = 0;
        ReverseFlag = 1;
        tau = 0;
    end
    
    %QN2 = zeros(((S1+N+1)+1)*(S2+N+1),((S1+N+1)+1)*(S2+N+1));
    %TmpV = zeros(((S1+N+1)+1)*(S2+N+1),1);
    %A{S2+N+1,S2+N+1} = [];
    %Modified_A{S2+N+1,S2+N+1} = [];
    
    %For the last A
    %TmpA = GetA(S2+N, S2+N, N, S1, S2, Lambda, mu1, mu2, tau, theta);
    %FinalDimension = length(find(sum(TmpA) ~= 1));
    
    for i = (S2+N): -1 : 0
        %Diagonal
        TmpA = GetA(i, i, N, S1, S2, Lambda, mu1, mu2, tau, theta, ReverseFlag);
        if i > 0
            TmpB = GetA(i - 1, i, N, S1, S2, Lambda, mu1, mu2, tau, theta, ReverseFlag);
        end
        if i < (S2 + N)
            TmpC = GetA(i + 1, i, N, S1, S2, Lambda, mu1, mu2, tau, theta, ReverseFlag);
        end
        
        if i == (S2+N)
            R = -TmpB/(TmpA);
        elseif i > 0
            load(['R_' num2str(i+1) '.mat']);
            R = -TmpB/(TmpA + R*TmpC);
        end
        save(['R_' num2str(i) '.mat'],'R');
    end
    A_00 = TmpA;
    A_10 = TmpC;
    load(['R_' num2str(1) '.mat']);
    ZeroConstraintMatrix = A_00 + R*A_10;
    %Ini Condition
    NormalizeConstraintVector = ones(size(R,1),1) + R*ones(size(R,2),1);
    R_Cum = R;
    for i = 2 : (S2+N)
        load(['R_' num2str(i) '.mat']);
        R_Cum = R_Cum*R;
        NormalizeConstraintVector = NormalizeConstraintVector + R_Cum*ones(size(R_Cum,2),1);
    end
    TmpNormalizeNum = max(NormalizeConstraintVector);
    ConstraintMatrix = [ZeroConstraintMatrix  NormalizeConstraintVector/TmpNormalizeNum];
    
    pi0 = ([zeros(1,size(ZeroConstraintMatrix,2)) 1/TmpNormalizeNum]) / full(ConstraintMatrix);
    
    piDimension = zeros(1,S2+N+1);
    piDimension(0+1) = length(pi0);
    piMatrix = zeros(S1+N+1,S2+N+1);
    piMatrix(:,1) = pi0;
    
    for i = 1 : S2+N
        load(['R_' num2str(i) '.mat']);
        prev_pi = piMatrix(1:piDimension((i-1)+1),(i-1)+1)';
        pi = prev_pi*R;
        piDimension(i+1) = length(pi);
        piMatrix(1:length(pi),i+1) = pi';
    end
    
    
    
    Distribution = piMatrix;
    if ReverseFlag == 1
        Distribution = Distribution';
    end
    %OldQN2 = QN2
    %save('Distribution_Ori.mat','Distribution')
    
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
    %PF: P {All servers in both stations are busy and the total number of people in the two queues = N}
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
    W1 = L1/(Lambda*p*(PE1 + PE12));
    W2 = L2/(Lambda*(1-p)*(PE2 + PE12));
    
    W_Test = p*W1 + (1-p)*W2;
    if isnan(W1) || isinf(W1)
        W_Test = W2;
    end
    if isnan(W2) || isinf(W2)
        W_Test= W1;
    end
    WaitingTime = W_Test-((1/mu1)*p+(1/mu2)*(1-p));
    %toc
    
    if ReverseFlag == 1
        %Reverse: Technical Reason
        %Reverse again because it is a global var.
        tmp = mu1;
        mu1 = mu2;
        mu2 = tmp;
        
    end
    return