InputFileName = ['input',num2str(1),'.txt'];
InitializeAndReadInput(InputFileName)
global Epsilon
global theta
global d1
global d2
global c1
global c2
global Budget
global Beta1
global Beta2
h1=[];
for i=1:7
    X=Acurr(i).x;
    tau=Acurr(i).fn;
    waiting=MeanWait_MatrixInverse(tau,X(1,1),X(1,2));
    max_securitylevel = d1*CalcR1(tau,theta) + d2*CalcR2(tau,theta);
    Cost = c1*X(1) + c2*X(2) + Beta1*CalcP(tau,theta) + Beta2*(1-CalcP(tau,theta));
    
    h1=[h1 ;[waiting-Epsilon,max_securitylevel,Cost-Budget,sqrt(Acurr(i).ConstraintCov)]];
end