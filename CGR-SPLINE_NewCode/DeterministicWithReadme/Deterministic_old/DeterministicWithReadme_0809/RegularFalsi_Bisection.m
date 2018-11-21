function [tau, meanWaitingTime,err, flag] = RegularFalsi_Bisection(lb, ub, S1, S2, Epsilon)
    %Note: Sol cannot < x0 (Due to the budget constraint)
    %flag:
    % 1 ==> Valid Solution
    %-1 ==> No Root
    Tol = 1e-7;
    MaxIter = 1000;
    
    f_lb = GetRemainTime(lb,S1,S2,Epsilon);
    f_ub = GetRemainTime(ub,S1,S2,Epsilon);
    CurSol = lb - ((ub - lb) / (f_ub - f_lb))*f_lb;
    f_CurSol = GetRemainTime(CurSol,S1,S2,Epsilon);
    for I = 2 : MaxIter
        if f_lb * f_CurSol < 0
            ub = CurSol;
            f_ub = f_CurSol;
        else
            lb = CurSol;
            f_lb = f_CurSol;
        end
        if (ub - lb)/((ub+lb)/2) < Tol
            tau = ub;
            meanWaitingTime = f_ub + Epsilon;
            err = (tau - lb)/lb;
            flag = 1;
            return
        end
        CurSol = lb - ((ub - lb) / (f_ub - f_lb))*f_lb;
        f_CurSol = GetRemainTime(CurSol,S1,S2,Epsilon);
        disp([ '[' num2str(CurSol) ',' num2str(f_CurSol) ']' ])
        if abs(CurSol - lb)/CurSol < Tol || abs(CurSol - ub)/CurSol < Tol
            %Use bisection if regular falsi failed
            CurSol = (lb + ub)/2;
            f_CurSol = GetRemainTime(CurSol,S1,S2,Epsilon);
        end
    end
    
    tau = ub;
    meanWaitingTime = f_ub + Epsilon;
    err = (tau - lb)/lb;
    flag = 1;
    return
end