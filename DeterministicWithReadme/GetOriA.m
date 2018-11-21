function A = GetOriA(i, j, N, S1, S2, Lambda, mu1, mu2, tau, theta)
    p = (-exp(-tau/theta) + 1)/(1 - exp(-1/theta));
    if i == j
        Dimension = N + S1 + 1 - max(0, i - S2);
        A = eye(N+S1+1);
        QN1 = zeros(Dimension,Dimension);
        for h1 = 0 : (Dimension-1)
            if h1 == 0
                QN1(h1+1,h1+1) = -p*Lambda;
                QN1(h1+1,(h1+1)+1) = p*Lambda;
            elseif h1 == (Dimension - 1)
                QN1(h1+1,h1+1) = -S1 * mu1;
                QN1(h1+1,(h1-1)+1) = S1 * mu1;
            else
                if h1 < S1
                    tmp = h1;
                else
                    tmp = S1;
                end
                QN1(h1+1,h1+1) = -(tmp * mu1 + p*Lambda);
                QN1(h1+1,(h1-1)+1) = tmp*mu1;
                QN1(h1+1,(h1+1)+1) = p*Lambda;
            end
        end
        UN = [eye(Dimension-1) zeros(Dimension-1,1);
              zeros(1,Dimension-1) 0];
        
        if i == N + S2
            Tmp = QN1 - min(i, S2)*mu2*eye(Dimension);
        else
            Tmp = QN1 - min(i, S2)*mu2*eye(Dimension) - (1-p)*Lambda*UN;
        end
        A(1:Dimension,1:Dimension) = Tmp;
        A = A(1:Dimension,1:Dimension);
        A = sparse(A);
        return
    elseif j > i %j = i+1
        if i < S2
            A = GetU(N,S1);
        else
            A = GetV(N-(i-S2),S1);
        end
        A = (1-p)*Lambda*A;
        A = sparse(A);
        return
    elseif j < i %j = i-1
        if i <= S2
            A = eye(N + S1 + 1);
        else
            A = GetV(N+1-(i-S2),S1)';
        end
        A = min(i,S2)*mu2*A;
        A = sparse(A);
        return
    end
    
    