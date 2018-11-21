function [flag, xstruct, iseed] = ThreeStageFlowline_c(problemparam, x, m, iseed)

[flagseed, perfmeas, perfmeasVar] = TSF(problemparam, x(1), x(2), x(3), x(4), m, iseed(1), iseed(2), iseed(3));
xstruct=struct('x', {}, 'fn', {}, 'FnVar', {}, 'FnGrad', {}, 'FnGradCov', {}, 'constraint', {}, 'ConstraintCov', {}, 'ConstraintGrad', {}, 'ConstraintGradCov', {});
xstruct(1).x             = x;
xstruct(1).fn            = perfmeas(1);
xstruct(1).constraint    = perfmeas(2) + TSFtrue(problemparam, [6;7;7;12]); %To make the rhs of inequality = 0
xstruct(1).FnVar         = perfmeasVar(1);
xstruct(1).ConstraintCov = perfmeasVar(2);

flag=flagseed(1);
iseed=[flagseed(2);flagseed(3);flagseed(4)];


end
