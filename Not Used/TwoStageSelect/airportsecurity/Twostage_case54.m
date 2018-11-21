global c1
global c2
global Beta1
global Beta2
global Budget
global Epsilon
c1 = 2;
c2 = 3;
Beta1 = 5;
Beta2 = 35;
Budget = 44;
Epsilon = 8;


global N
global lambda
global mu1
global mu2
global theta
N = 50;
lambda = 40;
mu1 = 8;
mu2 = 4;
theta = 0.0625;

global tautotal
global test
global average
global ba
global totalmobs
global totaltwo
global samplese
global big_s1
global big_s2
big_s1  = floor((Budget-Beta1)/c1);
big_s2 = floor((Budget-Beta1)/c2);
tau_matrix = xlsread('R&S_Case54_best_tau2.xls');
var_tau_matrix = xlsread('R&S_Case54_best_SE2.xls');
system_number = 133;
alpha = 0.5;
sample = 20 ;
[best_sys]=NSGS(var_tau_matrix,tau_matrix,sample,system_number,alpha);