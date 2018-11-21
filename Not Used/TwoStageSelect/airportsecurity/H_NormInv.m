function NormInv = H_NormInv(p)
    PLIM = 1E-18;
    P0 = -0.322232431088;
    P1 = -1;
    P2 = -0.342242088547;
    P3 = -0.0204231210245;
    P4 = -4.53642210148E-05;
    Q0 = 0.099348462606;
    Q1 = 0.588581570495;
    Q2 = 0.531103462366;
    Q3 = 0.10353775285;
    Q4 = 0.0038560700634;

    if p > 0.5
        p = 1 - p;
    end
    
    if p >= PLIM
        y = sqrt(-1 * log(p * p));
        vtemp = y + ((((y * P4 + P3) * y + P2) * y + P1) * y + P0) / ((((y * Q4 + Q3) * y + Q2) * y + Q1) * y + Q0);
    else
        vtemp = 8;
    end
        
    if p >= 0.5
        vtemp = -1 * vtemp;
    end
    
    NormInv = vtemp;