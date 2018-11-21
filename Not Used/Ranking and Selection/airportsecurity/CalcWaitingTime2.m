function WaitingTime = CalcWaitingTime2(S1,S2,tau)
    tic
    global OldQN2
    %N = 50;
    %S1 = 1;
    %S2 = 1;
    %Lambda = 36;
    %mu1 = 8;
    %mu2 = 6;
    %tau = 0.05;
    %theta = 0.0625;
    
    N = 50;
    %S1 = 3;
    %S2 = 4;
    Lambda = 18;
    mu1 = 8;
    mu2 = 6;
    %tau = 0.05;
    theta = 0.0625;
    
    
    %QN2 = zeros(((S1+N+1)+1)*(S2+N+1),((S1+N+1)+1)*(S2+N+1));
    %TmpV = zeros(((S1+N+1)+1)*(S2+N+1),1);
    %A{S2+N+1,S2+N+1} = [];
    %Modified_A{S2+N+1,S2+N+1} = [];
    
    %For the last A
    %TmpA = GetA(S2+N, S2+N, N, S1, S2, Lambda, mu1, mu2, tau, theta);
    %FinalDimension = length(find(sum(TmpA) ~= 1));
    
    
    for i = (S2+N): -1 : 0
        %Diagonal
        TmpA = GetOriA(i, i, N, S1, S2, Lambda, mu1, mu2, tau, theta);
        if i > 0
            TmpB = GetOriA(i - 1, i, N, S1, S2, Lambda, mu1, mu2, tau, theta);
        end
        if i < (S2 + N)
            TmpC = GetOriA(i + 1, i, N, S1, S2, Lambda, mu1, mu2, tau, theta);
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
    ConstraintMatrix = [ZeroConstraintMatrix  NormalizeConstraintVector];
    
    pi0 = ([zeros(1,size(ZeroConstraintMatrix,2)) 1]) / ConstraintMatrix;
    
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
    %OldQN2 = QN2
    save('Distribution_Ori.mat','Distribution')
    
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
    toc
    return