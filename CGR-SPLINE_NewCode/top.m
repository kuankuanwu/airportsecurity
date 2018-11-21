function top(problemseed, solverseed)
% Example: top([1;2;3],[2;3;4;5])
% Size of problemseed is specified by probparam(3)
% Size of solverseed  is the same as problem dimension (= 4 for the
% three-stage flowline problem)


xmin             = [1 1 1 1];
xmax             = [20 20 20 19];
budget           = 20000;
problemname      = 'ThreeStageFlowline_c'; 
                    % Replace with 'ThreeStageFlowline' if unable to 
                    % generate mex file for TSF.c. Note however that
                    % the oracle runs much slower with the .m file than 
                    % with the mex file. 

probparam        = [999 4 3 1 50 1000 20 20];
                    % probparam(1) = problem ID, not used
                    % probparam(2) = problem dimension = 4
                    % probparam(3) = nseeds = 3
                    % probparam(4) = nSecMeas = 1 
                    % probparam(5) = warm-up time
                    % probparam(6) = simulation end time
                    % probparam(7) = total service rate
                    % probparam(8) = total buffer space available
solvername       = 'cgRSPLINE_v2';
solverparam      = [1  50 1000 3.5 0.4 10 1.1 .002 0.1 1 8 1.1  500 1.5 0.95 0.65 0 5e-9];
                    % solverparam(1)= numfinalsols
                    % solverparam(2)= numrestarts
                    % solverparam(3)= kmax 
                    % solverparam(4)= q 
                    % solverparam(5)= delta 
                    % solverparam(6)= c1 
                    % solverparam(7)= c2
                    % solverparam(8)= cvthreshold
                    % solverparam(9)= biasthreshold
                    % solverparam(10)= yesgeom
                    % solverparam(11)= a in mk=a*c^k
                    % solverparam(12)= c in mk=a*c^k
                    % solverparam(13)= restart budget constant, a
                    % solverparam(14)= restart budget exponent, c (polynomial growth, a*r^c)
                    % solverparam(15)= alpha in (alpha_r = alpha*(1-c^(r+1)))
                    % solverparam(16)= c in (alpha_r = alpha*(1-c^(r+1)))
                    % solverparam(17)= yegeom (for b_r)
                    % solverparam(18)= ftol (function tolerance)


                    
%% Run cgR-SPLINE_v2 expnum number of times parallely
logfilename='TSF';
solverhandle=str2func(solvername);
[Awork, Abest, Acurr, Ainc, OPsolverseed] = solverhandle(xmin, xmax, ...
        problemname, probparam, problemseed, ...
        solverparam, solverseed, logfilename, ...
        budget); %#ok<ASGLU>

end
