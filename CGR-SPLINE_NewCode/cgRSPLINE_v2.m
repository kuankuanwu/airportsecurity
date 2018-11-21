%% ========================================================================
function [Awork, Acurr, Ainc, Abest, solverseed] = cgRSPLINE_v2(xmin, xmax, ...
    problemname, problemparam, problemseed, solverparam, solverseed, ...
    logfilename, budget)

global RefMatrixSize
global RefMatrix_tau
global LinearFalg
global mobs
mobs = 0;
RefMatrixSize = xmax;
global cal_tau
cal_tau=[];

%  ========================================================================    
%
%     Kalyani Nagaraj
%     Oklahoma State University
%     Last updated: February, 2017
%
%  ========================================================================    
%  INPUT
%       minX, maxX
%         range of X
%       problemname
%         problem name 
%       problemparam
%         vector of problem parameters 
%       problemseed
%         vector of input seeds to oracle
%       solverparam 
%          solverparam(1): numfinalsols
%          solverparam(2): numrestarts
%          solverparam(3): kmax 
%          solverparam(4): q 
%          solverparam(5): delta 
%          solverparam(6): c1 
%          solverparam(7): c2
%          solverparam(8): cvthreshold
%          solverparam(9): biasthreshold
%          solverparam(10): yesgeom
%          solverparam(11): constant, a
%          solverparam(12): base, c
%          solverparam(13): restart budget constant, a
%          solverparam(14): restart budget exponent, c (polynomial growth, a*r^c)
%          solverparam(15): alpha in (alpha_r = alpha*(1-c^(r+1)))
%          solverparam(16): c in (alpha_r = alpha*(1-c^(r+1)))
%          solverparam(17): yegeom (for b_r)
%       logfilename
%          solver log file
%       solverseed
%          siseed(1): seed to each restart loop
%       budg
%          total budget 
%
%  OUTPUT
%       Awork: a column vector with as many rows as there are restarts
%       Acurr: a data strucure (indexed by restart number) that
%              stores the current solution Y(r): 
%              Y(r).x, Y(r).fn, Y(r).FnVar, Y(r).constraint, 
%              Y(r).ConstraintCov, Y(r).lastRAipseed, Y(r).mk, Y(r).phat, 
%              Y(r).ghtrue, Y(r).FO
%       Ainc : data structure (indexed by restart number) that stores the
%              incumbent Z(r-1) (before comparision with the current soln.
%              Y(r) but at the updated sample size t_r):
%              Z(r-1).x, Z(r-1).fn, Z(r-1).FnVar, Z(r-1).constraint, 
%              Z(r-1).ConstraintCov, Z(r-1).lastRAipseed, Z(r-1).mk, Z(r-1).phat, 
%              Z(r-1).ghtrue, Z(r-1).FO
%       Abest: data structure (indexed by restart number) that stores the
%              incumbent Z(r) (after comparision with the current soln.
%              Y(r) ):
%              Z(r).x, Z(r).fn, Z(r).FnVar, Z(r).constraint, 
%              Z(r).ConstraintCov, Z(r).lastRAipseed, Z(r).mk, Z(r).phat, 
%              Z(r).ghtrue, Z(r).FO
%       solverseed: output solver seed
%  ========================================================================    
%  References:
%     [1] K. Nagaraj and R. Pasupathy,
%     Stochastically constrained simulation optimization on integer-ordered
%     spaces: the cgR-SPLINE algorithm, 
%     Under 2nd review with Operations Research, 2017
%
%     [2] H. Wang, R. Pasupathy, and B. Schmeiser,
%     Integer-Ordered Simulation Optimization using R-SPLINE: 
%     Retrospective Search with Piecewise-Linear Interpolation and
%     Neighborhood Enumeration, ACM TOMACS, 2013
%
%  ========================================================================
%     NOTES: 1. Introduced new variables 'repeat' and 'repcount' to track 
%               repeated restarts. This has been helpful with problems with
%               very small feasible regions relative to the search space, 
%               (that is the probability of generating an initial feasible 
%               solution very small, and GenInput invariably generates an 
%               infeasible solution.)

%            2. Code to check alpha_r feasibility of local solution X_r 
%               generated at the end of a restart is within cRSPLINE_v2.m
% 
%  ========================================================================


logfilename1=strcat(logfilename, '.txt');
logfid=fopen(logfilename1, 'w'); %General log file
logfilename2=strcat(logfilename, '_bestsols.txt');
logfid2=fopen(logfilename2, 'w'); %Log of the solutions returned at the end of each call to SPLINE for each of the restarts
logfilename3=strcat(logfilename, '_report.txt');
logfid3=fopen(logfilename3, 'w'); %The set M^*_r and the corresponding estimated obj func and constraint function values
logfilename4=strcat(logfilename, '_final.txt');
logfid4=fopen(logfilename4, 'w'); %Updates Z^*_r and returns its true g, h values 


% Problem parameters
id=int32(problemparam(2));
nseeds=int32(problemparam(3));
nsecMeas=int32(problemparam(4));
orchandle=str2func(problemname);

% Solver parameters
numrestarts=solverparam(2);  % size of A variables
yesgeom=solverparam(17);

% Declare o/p data structures
Awork=zeros(numrestarts, 1);

% Initialize data structures for Y_r, Z_r-1, Z_r:
% Solution from current call to R-SPLINE: Y_r
Acurr=struct('x',{},'fn',{},'FnVar',{},'constraint',{},'ConstraintCov',{},'lastRAipseed',{}, 'mk', {}, 'phat',{},'ghtrue',{},'FO',{});
%Incumbent Solution: Z_r-1
Ainc=struct('x',{},'fn',{},'FnVar',{},'constraint',{},'ConstraintCov',{},'lastRAipseed',{}, 'mk', {}, 'phat',{},'ghtrue',{},'FO',{});
%Current best after comparison with incumbent: Z_r
Abest=struct('x',{},'fn',{},'FnVar',{},'constraint',{},'ConstraintCov',{},'lastRAipseed',{}, 'mk', {}, 'phat',{},'ghtrue',{},'FO',{});

% Track total budget expended
totalcalls=0; 

  fprintf(logfid3, '\n--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------X_r-------------------------------\n');
    fprintf(logfid3, ' Restart                             Suggested                                                                         Current                                                              Total |                                Feasibility/\n');
    fprintf(logfid3, '     No.                              Solution                                                                            Best                                                               Work |                                  Optimality\n');
    fprintf(logfid3, '       r                                   X_r            gk             fk    phat     mk                               Z_r-1             gk             fk     phat   mkprev              (w_r) |             g             f\n');
    fprintf(logfid3, '--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------X_r-------------------------------\n');
    
  fprintf(logfid4, '\n-------------------------------------------------------------------------------------------------------------------|----------------------------------------------\n');
    fprintf(logfid4, ' Restart                              Suggested           g_kr           f_kr    phat       alpha_r     Total Work |             g              f     Feasibility/\n');
    fprintf(logfid4, '     No.                               Solution                                                                    |                                    Optimality\n');
    fprintf(logfid4, '     (r)                                  (Z_r)                                                              (w_r) |\n');
    fprintf(logfid4, '-------------------------------------------------------------------------------------------------------------------|----------------------------------------------\n');
global outter
% Initialize restart number r, counter for restarts    
starttime=clock;
r=1;
outter =1;
repcount=0;

Test_initial_point = [];
global total_sol
total_sol=[];
while(1)   % while work done <= budg
    fprintf(logfid, '\n\n\n\n======================== RESTART %d ===========================\n', r);
    fprintf(logfid2, '\n\n\n\n  RESTART %d\n', r);
    
    %% I. Generate a starting solution
    %[x0,solverseed]=GenInput(solverseed, xmin, xmax);% x0 is a column vector
    
    if LinearFalg ==0
       [x0,solverseed]=GenInput_AirportSecurity(solverseed, xmin, xmax);
    else 
       [x0,solverseed]=GetInput_Linear(solverseed,r,xmin,xmax,Acurr); 
    end 
    
    Test_initial_point(:,r) = x0;
    
    fprintf(logfid, '\nInitial solution = [');
    fprintf(logfid, '% 8d', x0);
    fprintf(logfid, ' ]\n\n');

    fprintf(logfid2, '  X_0 = [');
    fprintf(logfid2, '% 8d', x0);
    fprintf(logfid2, ' ]\n');

    %% II. Define b_r, alpha_r
    if yesgeom==0
        budgr=ceil(solverparam(13)*r^(solverparam(14))); %polynomial growth of restart budget
    else
        budgr=ceil(solverparam(13)*(solverparam(14)^r));   %geometric growth of restart budget
    end
    alphar=solverparam(15)*(1-solverparam(16)^(1+r));   %require X_r to be feas with prob at least alphar
    %alphar = solverparam(15)*(1-1/log(budgr));
    
    disp(['OuterIter:'   num2str(r) ',startpoint:' num2str(x0(1))   ',' num2str(x0(2)) ]);
    disp(['OuterIter:'   num2str(r) ',budgr:' num2str(budgr) ]);
    
    fprintf(logfid, '\nbr = % 8d, alphar = % 12.4f\n', budgr,alphar);
    fprintf(logfid2, '\nbr = % 8d, alphar = % 12.4f\n', budgr,alphar);

    
    %% III. Call cR-SPLINE
    repeat=0;  
    [ncalls, x1, x1phat, FO, dummy, mk, lastRAipseed] = cRSPLINE_v2(x0, ...
        orchandle, problemparam, problemseed, solverparam, ...
        logfid, logfid2, budgr, alphar);
    RefMatrix_tau = -1*ones(RefMatrixSize);
    
    problemseed = lastRAipseed;
    
    % k         : last RA iteraion number; needed to calculate m_{k_r}
    % lastiseed : input seed to the kth RA iteration;
    %             used to update ghat and fhat for backfill
    
    %% IV. Store and update X_r, Z_r-1, and Z_r.
    %Update Acurr for the current solution
    Acurr(r).x=x1.x;
    Acurr(r).fn=x1.fn;
    Acurr(r).FnVar=x1.FnVar;
    Acurr(r).constraint=x1.constraint;         %x1constraint is a col. vector
    Acurr(r).ConstraintCov=x1.ConstraintCov;   %x1ConstraintCov is a column vector
    Acurr(r).mk=mk;
    Acurr(r).phat=x1phat;
    Acurr(r).lastRAipseed=lastRAipseed;
    %Acurr(r).ghtrue=[x1.fn;TSFtrue(problemparam, x1.x)];  %****
    Acurr(r).FO=FO;
    
    totalcalls = totalcalls + ncalls;
    Awork(r)=totalcalls;    %cumulative work done so far
        
    if r==1    % The first restart
        % Update incumbent, Ainc.
        Ainc(r).x             = NaN*zeros(id, 1);
        Ainc(r).fn            = NaN;
        Ainc(r).FnVar         = NaN;
        Ainc(r).constraint    = NaN*zeros(nsecMeas, 1);
        Ainc(r).ConstraintCov = NaN*zeros(nsecMeas, 1);
        Ainc(r).lastRAipseed  = NaN*zeros(nseeds, 1);
        Ainc(r).mk            = NaN;
        Ainc(r).phat          = NaN; 
        Ainc(r).ghtrue        = NaN;
        Ainc(r).FO            ='NA'; % There isn't an incumbent
        
        
        Abest(r) = Acurr(r);

        if strcmp(Acurr(r).FO, '-')==0 % Current cR-SPLINE call returned a sample path feasible solution
            %prevlastiseed=Abest(r).lastRAipseed;
            % Not seeded since can simply access prev seed from
            % Abest(r).lastRAipseed
            % Print Report
            fprintf(logfid3, '% 8d   [', r);
            fprintf(logfid3, '% 8d', Acurr(r).x);
            fprintf(logfid3, ' ]% 14.6f% 15.9f% 8.4f% 7d% 103d |% 15s\n', Acurr(r).fn, Acurr(r).constraint, Acurr(r).phat(r), Acurr(r).mk, Awork(r), Acurr(r).FO);
        else                           
           fprintf(logfid3, '% 8d% 22s% 160d |% 43s\n', r, 'NONE', Awork(r), FO);
           repcount=repcount+1;
           repeat=1;
        end
        

    else  % If r>=2 
        
        % Update incumbent, Ainc.
        Ainc(r) = Abest(r-1);
        
        % POSSIBILITY I.
        % *** WILL NEVER OCCUR *** 
        % Previous RA iteration had to have returned a sample-path feasible solution
        % since restarts are repeated until a feasible solution is found.
        if strcmp(Acurr(r).FO, '-')==1 && strcmp(Ainc(r).FO, '-')==1 % no sample path feasible sol obtained so far
            Abest(r)=Ainc(r); % Doesn't matter what Abest is, since Acurr(r).FO='-'

            % Do nothing, simply log output 
            fprintf(logfid3, '% 8d% 22s% 131d |% 43s\n', r, 'NONE', Awork(r), Acurr(r).FO);

            repeat=1;          % Restart rth local seach from a different starting point 
            repcount=repcount+1;
        
            
            
        % POSSIBILITY II.
        % Previous restart returned a sample-path feasible solution, but the 
        % current restart did not. 
        elseif strcmp(Acurr(r).FO, '-')==1 && strcmp(Ainc(r).FO, '-')==0 
            Abest(r)=Ainc(r);

            % Print to report
            fprintf(logfid3, '% 8d                  NONE                                                               [', r);
            fprintf(logfid3, '% 8d', Abest(r).x);
            fprintf(logfid3, ' ]% 14.6f% 15.9f% 8.4f% 9d% 8d% 11d |% 43s\n', Abest(r).fn, Abest(r).constraint(1), Abest(r).phat, Abest(r).mk, 1, Awork(r), Acurr(r).FO);
            
            repeat=1;  %r=r-1;  % restart seach from a different starting point but at same sample size 
            repcount=repcount+1;

            
        % POSSIBILITY III.
        % *** WILL NEVER OCCUR ***    
        % First time a restart returns a sample-path feas. sol. (Will never
        % happen for r>=2 since a restart is repeated until a soln is
        % found.)
        elseif strcmp(Acurr(r).FO, '-')==0 && strcmp(Ainc(r).FO, '-')==1 

            Abest(r)=Acurr(r);
            
            fprintf(logfid3, '% 8d   [', r);
            fprintf(logfid3, '% 8d', Acurr(r).x);
            fprintf(logfid3, ' ]% 14.6f% 15.9f% 8.4f% 7d  % 20s% 65d |% 15s\n', Acurr(r).fn, Acurr(r).constraint(1), Acurr(r).phat, Acurr(r).mk, 'NONE', Awork(r), Acurr(r).FO);
        
            
        % POSSIBILITY IV.
        % Restart returns a sample-path feas. sol. 
        % Incumbent was also sample-path feasible.
        else
            
            if Acurr(r).mk>Ainc(r).mk % Update the incumbent
                
                fprintf(logfid, '\nCompare current solution with incumbent: Acurr(%d).mk=%d > Ainc(%d).mk=%d\n',r,Acurr(r).mk,r,Ainc(r).mk);
                fprintf(logfid, 'Call ORACLE at Ainc(%d).x with sample size %d and input seed Ainc(%d).lastRAipseed=%d,%d,%d,%d\n',r,Acurr(r).mk,r,Ainc(r).lastRAipseed);
                    
                % Call Oracle at Ainc(r).x with sample size Acurr(r).mk and
                % input seed = Ainc(r).lastRAipseed (not maintaining CRN)
                % Update other variables associated with Ainc(r).x
                fprintf(logfid, 'Size of Ainc(%d).x = [%d %d]\n',r, size(Ainc(r).x,1),size(Ainc(r).x,2));
                fprintf(logfid, 'Ainc(%d).x = [',r);
                fprintf(logfid, '%d ',Ainc(r).x);
                fprintf(logfid, ']\n\n');
                
                [flag, Ainc_r, ~] = orchandle(problemparam, Ainc(r).x, Acurr(r).mk, Ainc(r).lastRAipseed);
                Ainc(r).x=Ainc_r.x;
                Ainc(r).fn=Ainc_r.fn;
                Ainc(r).FnVar=Ainc_r.FnVar;
                Ainc(r).constraint=Ainc_r.constraint;
                Ainc(r).ConstraintCov=Ainc_r.ConstraintCov;

                if flag==0
                    totalcalls = totalcalls + mobs;
                    Awork(r)=totalcalls;
                end
                Ainc(r).mk=Acurr(r).mk;
                epsilon=min(max(sqrt(Ainc(r).ConstraintCov*Ainc(r).mk)/Ainc(r).mk^0.45, -Ainc(r).constraint), sqrt(Ainc(r).ConstraintCov*Ainc(r).mk)/Ainc(r).mk^0.1);
                delta=(log(sqrt(Ainc(r).ConstraintCov*Ainc(r).mk))-log(epsilon))/log(Ainc(r).mk); %delta is a vector of size nsecMeas
                delta(isnan(delta))=solverparam(5);
                if sum(Ainc(r).constraint > sqrt(Ainc(r).ConstraintCov*Ainc(r).mk)./(Ainc(r).mk).^delta)>0
                   Ainc(r).FO='-';        % Update Ainc(r).FO to read sample path infeasible 
                   
                end
                Ainc(r).phat=prod(tcdf((Ainc(r).mk).^(0.5-delta).*ones(nsecMeas,1)-Ainc(r).constraint./sqrt(Ainc(r).ConstraintCov), Ainc(r).mk-1)); 

                fprintf(logfid,'Update incumbent x=[ ' );
                fprintf(logfid, '%d ', Ainc(r).x);
                fprintf(logfid, '], new sample size = %d, Ainc(r).fn= %14.6f, Ainc(r).constaint= %.12f, Ainc(r).phat = %.12f, alpha_r=%.12f\n', Ainc(r).mk, Ainc(r).fn, Ainc(r).constraint(1), Ainc(r).phat, alphar);

            end

            if Acurr(r).mk<Ainc(r).mk % Update the current solution
                
                fprintf(logfid, 'Compare current solution with incumbent: Acurr(%d).mk=%d < Ainc(%d).mk=%d\n',r,Acurr(r).mk,r,Ainc(r).mk);
                fprintf(logfid, 'Call ORACLE at Acurr(%d).x with sample size %d and input seed Acurr(%d).lastRAipseed=%d\n',r,Ainc(r).mk,r,Acurr(r).lastRAipseed);
                    
                % Call Oracle at Acurr(r).x with sample size Ainc(r).mk and
                % input seed = Acurr(r).lastRAipseed (not maintaining CRN)
                % Update other variables associated with Acurr(r).x
                [flag, Acurr_r, ~] = orchandle(problemparam, Acurr(r).x, Ainc(r).mk, Acurr(r).lastRAipseed);
                Acurr(r).x=Acurr_r.x;
                Acurr(r).fn=Acurr_r.fn;
                Acurr(r).FnVar=Acurr_r.FnVar;
                Acurr(r).constraint=Acurr_r.constraint;
                Acurr(r).ConstraintCov=Acurr_r.ConstraintCov;
                
                if flag==0
                    totalcalls = totalcalls + mobs;
                    Awork(r)=totalcalls;
                end
                Acurr(r).mk=Ainc(r).mk;
                epsilon=min(max(sqrt(Acurr(r).ConstraintCov*Acurr(r).mk)/Acurr(r).mk^0.45, -Acurr(r).constraint), sqrt(Acurr(r).ConstraintCov*Acurr(r).mk)/Acurr(r).mk^0.1);
                delta=(log(sqrt(Acurr(r).ConstraintCov*Acurr(r).mk))-log(epsilon))/log(Acurr(r).mk); %delta is a vector of size nsecMeas
                delta(isnan(delta))=solverparam(5);
                if sum(Acurr(r).constraint > sqrt(Acurr(r).ConstraintCov*Acurr(r).mk)./(Acurr(r).mk).^delta)>0
                   Acurr(r).FO='-';        % Update Ainc(r).FO to read sample path infeasible 
                end
                Acurr(r).phat=prod(tcdf((Acurr(r).mk).^(0.5-delta).*ones(nsecMeas,1)-Acurr(r).constraint./sqrt(Acurr(r).ConstraintCov), Acurr(r).mk-1)); 

                fprintf(logfid,'Update current x=[ ' );
                fprintf(logfid, '%d ', Acurr(r).x);
                fprintf(logfid, '], new sample size = %d, Acurr(r).fn=%14.6f, Acurr(r).constaint= %.12f, Acurr(r).phat = %.12f, alpha_r=%.12f\n', Acurr(r).mk, Acurr(r).fn, Acurr(r).constraint(1), Acurr(r).phat, alphar);
            end

            % Print Report
            if strcmp(Ainc(r).FO,'-')==1 % If incumbent gets updated as sample-path infeasible
                inc_feas=0;
            else
                inc_feas=1;
            end
            fprintf(logfid3, '% 8d   [', r);
            fprintf(logfid3, '% 8d', Acurr(r).x);
            fprintf(logfid3, ' ]% 14.6f% 15.9f% 8.4f% 7d   [', Acurr(r).fn, Acurr(r).constraint(1), Acurr(r).phat, Acurr(r).mk);
            fprintf(logfid3, '% 8d', Ainc(r).x);
            fprintf(logfid3, ' ]% 14.6f% 15.9f% 8.4f% 9d% 8d% 11d |% 15s\n', Ainc(r).fn, Ainc(r).constraint(1), Ainc(r).phat, Ainc(r).mk, inc_feas, Awork(r), Acurr(r).FO);

            % Update Abest:
            if strcmp(Ainc(r).FO,'-')==1 || Ainc(r).phat < alphar % After updating incumbent, 
                                                                  % if incumbent is sample-path infeasible 
                                                                  % or not alpha_r feasible
                inc_good=0;
            else
                inc_good=1;
            end
            
            if strcmp(Acurr(r).FO,'-')==1 || Acurr(r).phat < alphar % After updating current soln, 
                                                                    % if it is sample-path infeasible 
                                                                    % or not alpha_r feasible
                curr_good=0;
            else
                curr_good=1;
            end
            % Note that inc_good and curr_good can never be 0 at the same
            % time since only one of the two is ever updated. The other is
            % always sample-path as well as alphar feasible
            
            
            if inc_good == 1 && curr_good == 0
                Abest(r)=Ainc(r);
            elseif inc_good == 0 && curr_good == 1
                Abest(r)=Acurr(r);
            elseif inc_good == 0 && curr_good == 0 % Can never be possible!
                Abest(r)=Ainc(r);                  % Error: set as default
            else
                if Acurr(r).fn<Ainc(r).fn
                    Abest(r)=Acurr(r);
                else
                    Abest(r)=Ainc(r);
                end
            end           
            
            
        end
        
    end  
    
    % Print Final Report
    if strcmp(Abest(r).FO, '-') == 0 % R-SPLINE found a sample path feasible solution
        % Print Report
        fprintf(logfid4, '% 8d     [', r);
        fprintf(logfid4, '% 8d', Abest(r).x);
        fprintf(logfid4, ' ]% 14.6f% 15.9f% 8.4f% 14.4f% 15d |% 16s\n', Abest(r).fn, Abest(r).constraint(1), Abest(r).phat, alphar, Awork(r), Abest(r).FO);
        % End report
    else
        fprintf(logfid4, '% 8d                 NONE                                                        % 14.4f% 15d |% 45s\n', r, alphar, Awork(r), Abest(r).FO);
    end

    %% V. End restarts if total budget exceeds bmax
    if totalcalls>budget
       break;
    end
    
    % Increment r if FO ~= '-' (that is if restart returns a sample-path feasible solution), else restart rth local search with a new seed 
    if repeat==0 
        r=r+1;
        outter=outter+1;
        repcount=0;
    end
    
    
end
endtime=clock;
fprintf(logfid, 'Total elapsed time (in minutes) = %0.12f\n\n', double(etime(endtime, starttime))/60);

save(strcat(logfilename, '_vars.mat'), 'Acurr', 'Ainc', 'Abest', 'Awork', 'Test_initial_point','total_sol','cal_tau');
    
fclose(logfid);
fclose(logfid2);
fclose(logfid3);
fclose(logfid4);
    
end




%%  =======================================================================
function [x0, seed] = GenInput(seed, xmin, xmax)

%   NOTE: x0 is a column vector
%   Written for the three-stage flowline problem.
%   Requires as many seeds as the dimension of the search space.
%   =======================================================================
id=length(xmin);
x0=zeros(id,1);
for i=1:id
	[seed(i), u] = u16807d(seed(i));
    x0(i) = floor(xmin(i) + u*(xmax(i) - xmin(i)) + 0.5);
end
end



%% ========================================================================
function [iseed,u]=u16807d(iseed)
% ========================================================================
%     bruce schmeiser     january 1989.                                   .
%     a linear congruential pseudorandom number generator                 .
%       using constant 16807 and modulus (2**31)-1.                       .
%              iseed = iseed*16807 (mod 2^31 -1)                          .
%     for implementations that don't require double precision, see        .
%       s.k. park and k.w. miller, "random numbers generators: good       .
%         ones are hard to find," cacm, 31, 10 (October 1988), 1192-1201. .
%     in correct implementations, starting with a seed value of 1 will    .
%       result in a seed value of 1043618065 on call number 10001.        .
%..........................................................................
%     input:  iseed.   integer.                                           . 
%                        chosen from [1,2147483646] on the first call.    .
%                        thereafter, the value returned from the last call.
%     output: iseed.   integer.                                           .
%                        to be used in the next call.                     .
%     output: u16807d.  real.                                             .
%                        a pseudorandom number in (0,1).                  .
%..........................................................................
u=0;
while u<=0 || u>=1
    iseed = mod (iseed * 16807,2147483647);
    u = iseed / 2147483648;
end
end

function [x0, seed] = GenInput_AirportSecurity(seed, xmin, xmax)
%Bigghost Modified
global Budget
global Beta1
global Beta2
global c1
global c2


id=length(xmin);
x0=zeros(id,1);
[u, seed] = mrg32k3a(seed);
x0(1) = floor(xmin(1) + u*(xmax(1) - xmin(1)) + 0.5);
Upper = floor((Budget - Beta1 - c1*x0(1))/c2);
Lower = floor(max((Budget - Beta2 - c1*x0(1))/c2,0));
% [seed, u] = u16807d(seed);
[u, seed] = mrg32k3a(seed);
x0(2) = floor(Lower + u*(Upper - Lower) + 0.5);

% for i=1:id
% 	[seed, u] = u16807d(seed);
%     x0(i) = floor(xmin(i) + u*(xmax(i) - xmin(i)) + 0.5);
% end
end

function [x0,seed] = GetInput_Linear(seed,r,xmin,xmax,Acurr)
  
global Budget
global Beta1
global Beta2
global c1
global c2

 id=length(xmin);
 x0=zeros(id,1);

 td=0;

if r < 3
    
[seed,u] = u16807d(seed);
% [u, seed] = mrg32k3a(seed);
x0(1) = floor(xmin(1) + u*(xmax(1) - xmin(1)) + 0.5);
Upper = floor((Budget - Beta1 - c1*x0(1))/c2);
Lower = floor(max((Budget - Beta2 - c1*x0(1))/c2,0));
[seed,u] = u16807d(seed);
% [u, seed] = mrg32k3a(seed);
x0(2) = floor(Lower + u*(Upper - Lower) + 0.5);
 
end   


if r > 2 
    
  if   Acurr(r-1).x(1) == Acurr(r-2).x(1) && Acurr(r-1).x(2) == Acurr(r-2).x(2)
     [seed,u] = u16807d(seed);
     x0(1) = floor(xmin(1) + u*(xmax(1) - xmin(1)) + 0.5);
     Upper = floor((Budget - Beta1 - c1*x0(1))/c2);
     Lower = floor(max((Budget - Beta2 - c1*x0(1))/c2,0));
     [seed,u] = u16807d(seed);
     x0(2) = floor(Lower + u*(Upper - Lower) + 0.5);
      
  else
  td = (Acurr(r-1).x(2)-Acurr(r-2).x(2))/ (Acurr(r-1).x(1)-Acurr(r-2).x(1)); 
  
  [seed,u] = u16807d(seed);
  x0(1) = floor(xmin(1) + u*(xmax(1) - xmin(1)) + 0.5);
 
  
  x0(2)= floor(Acurr(r-1).x(2)+ td*(x0(1)-Acurr(r-1).x(1))+0.5);
     Upper = floor(max((Budget - Beta1 - c1*x0(1))/c2,0));
     Lower = floor(max((Budget - Beta2 - c1*x0(1))/c2,0));
  if x0(2)< Lower 
     x0(2) =  Lower;
  end
  if x0(2)>  Upper 
     x0(2) = Upper;
  end
  end    
end
end

