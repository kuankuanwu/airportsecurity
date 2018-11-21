function TMP = H_ChiPdf(N, C, lngam) 
    FLN2 =  N/2;
    TMP = -FLN2 * log(2) - lngam(N) + (FLN2-1)*log(C) - C/2;
    TMP = exp(TMP);
end