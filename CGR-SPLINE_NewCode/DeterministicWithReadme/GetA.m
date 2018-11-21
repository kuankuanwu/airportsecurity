function A = GetA(i, j, N, S1, S2, Lambda, mu1, mu2, tau, theta)
%  purpose: To evaluate the transition rate matrix in non-selectee lane given the 
%           current state and the future state of selectee lane.
%           The output is always a (N+S1+1)*(N+S1+1) matrix, not in dynamic
%           size. If the dimension of the evaluated matrix is less than
%           N+S1+1, it will fill the identical matrix at the right-bottom
%           side of the evaluated matrix as the output.
%  input
%   i:          The current state of the selectee lane. (i.e. The current
%               queue length plus the number of visitors
%               being processed in selectee lane)
%   j:          The future state of the selectee lane. (i.e. The future
%               queue length plus the number of visitors
%               being processed in selectee lane)
%   N:          The maximum queue length limit of the sumation of queue
%               lengths in non-selectee lane and selectee lane.
%   S1:         The number of servers in non-selectee lane.
%   S2:         The number of servers in selectee lane.
%   Lambda:     The arrival rate of visitors. (visitors/hr)
%   mu1:        The service rate per server in non-selectee lane
%               (visitors/hr)
%   mu2:        The service rate per server in selectee lane
%               (visitors/hr)
%   tau:        The threshold to decide the visitors go to non-selectee lane
%               or selectee lane. If the risk value of the visitor <= tau, this
%               visitor would be assigned to non-selectee lane otherwise to selectee lane.
%   theta:      Parameters of truncaed exponential distribution which is
%               assumed as the distriubtion of visitors' risk value
%  output
%    A:         To evaluate the transition rate matrix in non-selectee lane given the 
%               current state and the future state of selectee lane.
%  variables:
%   mu1:        The service rate per server in non-selectee lane
%               (visitors/hr)
%   mu2:        The service rate per server in selectee lane
%               (visitors/hr)
%   theta:      Parameters of truncaed exponential distribution which is
%               assumed as the distriubtion of visitors' risk value

    p = (-exp(-tau/theta) + 1)/(1 - exp(-1/theta));
    if i == j
        %No change in selectee lane
        Dimension = N + S1 + 1 - max(0, i - S2);
        A = eye(N+S1+1);
        QN1 = zeros(Dimension,Dimension);
        if Dimension > 1
            %Regular: The non-selectee lane exists
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
        else
            %Special case: The non-selectee lane does not exist
            QN1 = 0;
        end
        UN = [eye(Dimension-1) zeros(Dimension-1,1);
              zeros(1,Dimension-1) 0];
        
        if i == N + S2
            %Margin: Current non-selectee lane is full (i.e. No visitors can
            %enter non-selectee lane during this time)
            Tmp = QN1 - min(i, S2)*mu2*eye(Dimension);
        else
            if i < S2
                %If idle server in selectee lane exists, system allows visitors
                %which risk value greater than tau go to selectee lane
                %directly even the non-selectee lane is full
                Tmp = QN1 - min(i, S2)*mu2*eye(Dimension) - (1-p)*Lambda*eye(Dimension);
            else
                %If idle server in selectee lane does not exist, system does not allows visitors
                %which risk value greater than tau go to selectee lane
                %if the non-selectee lane is full
                Tmp = QN1 - min(i, S2)*mu2*eye(Dimension) - (1-p)*Lambda*UN;
            end
        end
        A(1:Dimension,1:Dimension) = Tmp;
        return
    elseif j > i %j = i+1
        %The queue length and the number of visitors in the selectee lane
        %increase by 1
        A = zeros(N + S1 + 1, N + S1 + 1);
        if i < S2
            Dimension = S1 + N + 1;
            for h1 = 0 : (Dimension-1)
                A(h1+1,h1+1) = 1;
            end
        else
            Dimension = S1 + N + 1 - (i - S2);
            for h1 = 0 : (Dimension-1-1)
                A(h1+1,h1+1) = 1;
            end
        end
        
        A = (1-p)*Lambda*A;
        return
    elseif j < i %j = i-1
        %The queue length and the number of visitors in the selectee lane
        %decrease by 1
        if i <= S2
            A = eye(N + S1 + 1);
        else
            A = zeros(N + S1 + 1, N + S1 + 1);
            Dimension = N + S1 + 1 + 1 - (i-S2);
            for h1 = 0 : (Dimension-1-1)
                A(h1+1,h1+1) = 1;
            end
        end
        A = min(i,S2)*mu2*A;
        return
    end
    
    