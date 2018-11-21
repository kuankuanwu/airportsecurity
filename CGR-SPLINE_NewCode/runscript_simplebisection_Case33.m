%%============TOP LEVEL CODE FOR CALLING SOLVER CODE ======================
%
% *** DATE
%        Last updated August, 2014
%
% *** AUTHOR
%        Kalyani Nagaraj
%        Virginia Tech, Blacksburg, VA 
%        kalyanin AT vt DOT edu
%
% *** NOTE
%        The following is an example configuration+script file for solving  
%        the stochastically constrained three-stage flowline problem using 
%        cgR-SPLINE. The user is required to change entries under "Problem
%        Parameters" to intialize a different problem.
%
%==========================================================================
%  INPUT
%        x0 
%              Matrix of size = 'numinitsols' X 'dim' of 
%              'numinitsols' initial solutions to the solver.
%              Each initial solution is of size 1 X 'dim'.
%              Num rows of x0 determined by solver requirements. 
%              Num cols of x0 determined by problem instance (for
%              example, 10-bus scheduling problem has 'dim'=10 whereas the
%              20-bus scheduling problem has 'dim'=20). 
%              variable type, that is integer/real/mixed, determined by 
%              problem type.
%        budget
%              Matrix of size 1 X 'numfinalsols' where 'numfinalsols' is 
%              the number of solutions (or sets of solutions) returned by 
%              the solver.
%              For example, if budget = [500 1000] then the solver is 
%              expected to return a total of two solutions (or sets of 
%              solutions), once every 500 calls to the oracle.
%              Whether num cols equals 1 or can exceed 1 is determined 
%              by the solver's capabilities.
%              Integer valued.
%        problemname
%              Problem function name.
%              Character string.
%        problemparam
%              A vector of problem specific parameters.
%              probparam = ['binary', 'dim', nprimMeas, nsecMeas, <...>]
%              'binary'   : yes=1/no=0 
%              'dim'      : Dimension of decision variable.
%                           For example, if problem is in R^d, 'dim'=d.
%                           Integer valued.
%              'nprimMeas': Number of objective functions.
%              'nsecMeas' : Number of secondary functions.
%              <...>      : A vector of other problem specific parameters.
%                           For example, for the integer-ordered (s,S) 
%                           Inventory Problem this list includes 
%                           'fixedordercost', 'perunitcost', 'holdingcost', 
%                           'warmupperiodlength'.
%                           Integer/real/mixed depending on problem.
%        solvername
%              Solver function name.
%              Character string.
%        solverparam
%              solverparam =['numinitsols', 'numfinalsols', <...>]
%              'numinitsols' : Number of initial solutions to solver
%              'numfinalsols': Number of solutions (or sets of solutions)
%                              returned by solver. 
%              <...>         : A vector of other solver specific paramaters.
%                              For example, for cR-SPLINE: 'kmax', 'q', 'm0',
%                              'b0', 'c1', 'c2'. 
%                              Integer/real/mixed depending on solver.
%        problemseed
%              A vector of size 1 X 'nprobseeds' of starting seeds or  
%              substream indices
%        problemseedparam
%              problemseedparam = ['nprobseeds', 'pstreamindex', <...>]
%              'nprobseeds'   : Number of seeds/substreams required by
%                               problem.
%                               Integer valued.
%              'pstreamindex' : 1, if problemseed is a vector of substream
%                               indices, and
%                               0, if problemseed is a vector of seeds.
%              <...>          : Other parameters
%                               Integer/real/mixed.
%        solverseed
%              vector of size 1 X 'nsolvseeds' of starting seeds or substream
%              indices
%        solverseedparam
%              solverseedparam = ['nsolvseeds', 'sstreamindex', <...>]
%              'nsolvseeds'  : Number of seeds/substreams required by
%                              solver.
%                              Integer valued
%              'sstreamindex': 1, if solverseed is a vector of substream
%                              indices, and
%                              0, if solverseed is a vector of seeds.
%              <...>         : Other paramaters
%                              Integer/real/mixed.
%        logfilename
%              A log file base name.
%
% OUTPUT
%        ncalls
%              Total number of calls made to the oracle. 
%        x
%              An array of size = 'numfinalsols' X 'dim' X 'maxsolsetsize' 
%              of solutions returned by solver.
%              Note that integer variable 'maxsolsetsize' is generated by 
%              the solver. It denotes the maximum number of estimated best 
%              solutions returned after budget(1), budget(2), ..., 
%              budget(numfinalsols) number of calls to the oracle.
%        Fn 
%              A matrix of size 'numfinalsols' X ('nprimMeas'+'nsecMeas') X
%              'maxsolsetsize' of performance measure estimates.
%              fn = [g1, g2, ..., g_nprimMeas, h1, h2, ..., h_nsecMeas]
%              g1, g2, ..., g_nprimMeas are primary (or objective) function 
%              estimates. 
%              h1, h2, ..., h_nsecMeas are secondary (or constraint)
%              function estimates. 
%              Reals.
%        FnVar
%              Array of covariance estimators of Fn.
%              Reals. 
%        Grad     
%              Gradient estimate of the true perfomance measures at x
%              Will typically be a matrix of size 
%              ('nprimMeas'+'nsecMeas') X 'dim',
%              unless 'numfinalsols' > 1, or maxsolsetsize > 1, or both.
%              Reals.
%        GradVar
%              Array of covariance estimators of Grad. 
%              Reals.
%        problemseed     
%              Output problem seed(s).
%        solverseed     
%              Output solver seed(s).
% =========================================================================


%% SET PROBLEM AND SOLVER PARAMETERS:
addpath('problemcode');
addpath('airportsecurity');  
%% Problem Paramaters
global log_test
%(Constraint) % (Epsilon lambda mu1 mu2 are the same unit)
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
Budget = 95;
Epsilon = 15;

%(Waiting Time)  
global N
global lambda
global mu1
global mu2
global theta
N = 50;
lambda = 40/60;
mu1 = 10/60;
mu2 = 1/60;
theta = 0.55;

m0 = 2000;
geomc = 3;

global AirportSecurityFlag
AirportSecurityFlag = 1;
global StartificationFlag
StartificationFlag = 0;
global LinearFalg
LinearFalg = 1;
global TSFtrueFlag
TSFtrueFlag = 0;

for trial = 2:2
log_test_name =['Debug_log_SimpleBisection_Case33_smallbase2_RandomInitial_trial',int2str(trial),'.txt']; 
log_test=fopen(log_test_name, 'w'); 
problemname      = 'AirportSecurity';
problemparam     = [0 2 4 1 200 0.001 0.005 0];
                             % [binary dim nprimMeas nsecMeas ...]
                             % see ThreeStageFlow.m for details
problemseedparam = [4, 0] ;  % nprobseeds = 3
                             % pstreamindex = 0 implies ThreeStageFlow.m 
                             % requires actual seeds values
problemseed      = [11111111,22222222,33333333,44444444];  
xmin             = [0 0];
xmax             = [ceil((Budget-Beta1)/c1) ceil(Budget/c2)];  
                             % [min; max]
                             % cgRSPLINE does not require an initial sol 
                             % but the interval of search space X 
                             
%% Solver Paramaters

solvername       = 'cgRSPLINE_v2';
solverparam      = [1 100 1000 3.5 0.4 8 10 .002 0.1 1 m0 geomc ...
                    500000 1.1 0.05 0.65 1 0]; 
%solverparam      = [1 50 1000 3.5 0.4 8 10 .002 0.1 1 8 1.1 ...
%                    500 1.5 0.95 0.65];  
                             % [numfinalsols numrestarts kmax q delta ... ]
                             % see cgRSPLINE.m for details
solverseedparam  = [1, 0];   % nsolvseeds = 1
                             % sstreamindex = 0 implies cgRSPLINE.m 
                             % requires actual seeds values
solverseed       = 872914+541311*trial;    %48906,10000,20000
logfilename      = strcat('cgRSPLINE-on-AirportSecurity_SimpleBisection_Case33_new_Linear_',int2str(solverseed),'small_base2');  
                             % log file base name 
budget           = 100000000;    % 'numfinalsols'=1


                                                         
%% CALL cgR-SPLINE
solverhandle=str2func(solvername);
[Awork, Acurr, Ainc, Abest, solverseed] = solverhandle(xmin, xmax, ...
    problemname, problemparam, problemseed, solverparam, solverseed, ...
    logfilename, budget);

end