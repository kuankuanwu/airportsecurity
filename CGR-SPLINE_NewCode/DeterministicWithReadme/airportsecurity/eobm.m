function vybarhat = eobm(n,y,m)
%  bruce schmeiser.  august 4, 2004.
%  purpose: to compute the obm estimate to estimate var(ybar)
%  reference:  chapter ii.4 in wheyming tina song's dissertation
%  parameter definitions
%  input
%    n:     number of observations
%    y:     vector of observations
%    m:     batch size, an algorithm parameter (0 .lt. m .lt. n)
%  output
%    vybarhat: estimated variance of the sample mean
% 
%  variable definitions
%    i:     index for both data points and batches
%    sumy:  sum of data points 1 to i
%    sumd:  sum of data points 1 to i-m
%    bsum:  sum of data points i-m+1 to i   (batch sum)
%    sum:   sum of the batch sums
%    sum2:  sum of the squared batch sums
%    ybar:  sample average
%    sumbm: sum of the batch means
%    sum2bm:sum of the squared batch means

    %process the first m observations---the first complete batch        
    sumy = 0;
    for i = 1 : m
        sumy = sumy + y(i);
    end
    sum = sumy;
    sum2 = sum*sum;
    %process observations m+1 through n
    sumd = 0;
    for i = (m+1) : n
        sumd = sumd + y(i-m);
        sumy = sumy + y(i);
        bsum = sumy - sumd;
        sum = sum + bsum;
        sum2 = sum2 + bsum*bsum;
    end
    %convert from batch sums to batch means
    ybar = sumy / n;
    sumbm = sum / m;
    sum2bm = (sum2/m) / m;
    %calculate the obm estimator
    vybarhat = ((sum2bm - ybar * (2*sumbm - (n-m+1)*ybar))) / ...
        (((n-m+1)*(n-m))/m);
    return
end