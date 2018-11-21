function [expo, iseed] = rexpo(mean,iseed)
    %purpose: To generate an observation of an exponential distribution
    %with expected value equal to mean
    %input: mean (expected value of the exponential distribution)
    [u, iseed] = mrg32k3a(iseed);
    expo = -log(1 - u)*mean;
end