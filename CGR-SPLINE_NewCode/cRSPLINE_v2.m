%% ========================================================================

function [ncalls, x1, x1phat, FO, iseed, mk, lastRAipseed] = ...
    cRSPLINE_v2(x0, orchandle, problemparam, iseed, solverparam, ...
    logfid, logfid2, budget, alphar)
    global cal_tau
    global trial
%  ========================================================================
%
%  Kalyani Nagaraj
%  Oklahoma State University
%  Last updated: February, 2017
%
%  ========================================================================    
%  INPUT   (self explanatory)
%
%  OUTPUT  (self explanatory)
%
%  ========================================================================    
%  REFERENCES:
%
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
%  NOTES: 1. cR-SPLINE_v2 does not prematurely terminate restarts.
%         2. cR-SPLINE_v2 guarantees by alpha_r feasibility of final
%            soln. Y_r by tracing back the trajectory. 
% 
%  ========================================================================

% Problem parameters 
% id=int32(problemparam(2));
% nseeds=int32(problemparam(3));
% nsecMeas=int32(problemparam(4));

% Solver parameters
kmax=solverparam(3);
q=solverparam(4);
%delta=solverparam(5)*ones(nsecMeas,1);
c1=solverparam(6);
c2=solverparam(7);
% cvthresh=solverparam(8);
% biasthresh=solverparam(9);
yesgeom=solverparam(10);
geoma=solverparam(11);
geomc=solverparam(12);
%ftol=solverparam(18);

global RefMatrix_tau
global RefMatrix_AvgWT
global RefMatrix_Se_AvgWT
global RefMatrix_SecurityLevel
global RefMatrix_Cost
global RefMatrix_mobs
global RefMatrixSize
global total_sol 
global TSFtrueFlag
global RefMatrix_Seed
% Declare global variables
global trajectory;
trajectory=[]; % Initialize empty trajectory at the start of each restart
global k

% Initialize cR-SPLINE variables
x1=struct('x', {}, 'fn', {}, 'FnVar', {}, 'FnGrad', {}, 'FnGradCov', {}, 'constraint', {}, 'ConstraintCov', {}, 'ConstraintGrad', {}, 'ConstraintGradCov', {});
ncalls=0;	% tracks the total calls made to the oracle
x1(1).x=x0;
k=1;


% PRINT TO LOG2
fprintf(logfid2, '-------------------------------------------------------------------------------------------------------------------------|--------------------------------------------\n');
fprintf(logfid2, '  k    delta     m_k       w_k                                 X_k           g_k  se(g_k)         f_k             se(f_k)|             g           f       Feas    Opt\n');
fprintf(logfid2, '-------------------------------------------------------------------------------------------------------------------------|--------------------------------------------\n');
% END PRINT TO LOG2
        
% countsame=0;  %For premature termination of restarts
%. Begin Retrospective Iterations
disp(['Start RSPLINE: ((S1,S2), iseed, ncalls) = ' num2str(x0(1)) ',' num2str(x0(2)) ',' num2str(iseed(1)) ',' num2str(iseed(2)) ',' num2str(iseed(3)) ',' num2str(ncalls)])
while (1)  % while ncalls<budget && k<kmax
    fprintf(logfid, '==== RETROSPECTIVE ITERATION k = %d ====\n', k);
  
    % Calculate mk, bk.
    if yesgeom==1
        mk=ceil(geoma*geomc^(k-1));
        bk=ceil(c1*c2^k);
    else
        mk=c1*ceil(k^q);
        bk=c2*ceil(k^q);
    end
    
    fprintf(logfid, 'mk = %d, bk = %d\n', mk, bk);
    fprintf(logfid, '===== BEGIN SPLINE =====\n');
    numbers=0;
    xk=x1.x;  % Initial sol. = sol. retured by previous call to SPLINE 
    lastRAipseed=iseed;  % Maintain a copy of i/p seed to the ongoing inner iteration 
    
    disp(['Call SPLINE: New Sample Path'])
    disp(['InnerIter:' num2str(k) ',mk=' num2str(mk)])
    disp(['InnerIter:' num2str(k) ',bk=' num2str(bk)])
    disp(['InnerIter:' num2str(k) ',seed:' num2str(iseed(1)) ',' num2str(iseed(2))  ',' num2str(iseed(3))  ])
    
    RefMatrix_tau = -1*ones(RefMatrixSize);
    RefMatrix_AvgWT = -1*ones(RefMatrixSize);
    RefMatrix_Se_AvgWT = -1*ones(RefMatrixSize);
    RefMatrix_SecurityLevel = -1*ones(RefMatrixSize);
    RefMatrix_Cost = -1*ones(RefMatrixSize);
    RefMatrix_mobs = -1*ones(RefMatrixSize);
    RefMatrix_Seed = -1*ones(1,4,RefMatrixSize(1),RefMatrixSize(2));
    Solutions=zeros(1,kmax);
    [splinencalls, stackflag, x1, iseed] = SPLINE(orchandle, ...
        problemparam, solverparam, xk, mk, bk, k, iseed, logfid);	
    ncalls = ncalls + splinencalls;
    
    %calculate the solutions number 
    if stackflag==0
        numbers = numbers + splinencalls/(mk+problemparam(5));
        Solutions(k)=numbers;
        disp(['numbers of solution =' num2str(numbers) 'in' num2str(k)])
    end 
    total_sol =[total_sol;Solutions];
    % Calculate phat=P{X^*\in \mathcal{F}}. Formula of phat assumes independent constraint functions
    epsilon=min(max(sqrt(x1.ConstraintCov*mk)/mk^0.45, -x1.constraint), sqrt(x1.ConstraintCov*mk)/mk^0.1);
    delta=(log(sqrt(x1.ConstraintCov*mk))-log(epsilon))/log(mk); %delta is a vector of size nsecMeas
    delta(isnan(delta))=solverparam(5);    
    x1phat=prod(tcdf(mk.^(0.5-delta)-x1.constraint./sqrt(x1.ConstraintCov), mk-1));
    
    disp(['InnerIter:' num2str(k) 'end'])
    disp(['InnerIter:' num2str(k) ',endpoint:'  num2str(x1.x(1)) ',' num2str(x1.x(2)) ])
    disp(['InnerIter:' num2str(k) ',ncalls' num2str(ncalls) ])

    FO='-';                         % FO defaults to 'sample-path infeasible' 
    if stackflag == 1               % stackflag=1 implies no feasible solution found.
        if k==1
            fprintf(logfid2, 'Initial solution is infeasible!!!\n');
        else
            fprintf(logfid2, 'All solutions in stack are infeasible!!!\n');
        end
        return                      % Return control to cgR-SPLINE, restart local search
    end  
    fprintf(logfid,'\nSPLINE returned x with phat=%.12f when alphar=%.12f\n',x1phat,alphar);
    fprintf(logfid,'\n===== SPLINE ENDED =====\n');
   
    % Trace sample-path backwards for an alpha_r-feasible solution      
    if (k>=kmax || ncalls>=budget) && x1phat < alphar  % last inner loop and alphar not attained
                                                     % first part checks
                                                     % termination criteria
                                                     % second part evaluates alphat
        fprintf(logfid, '\n\n===== BEGIN BACKTRACK =====\n');                                           
        fprintf(logfid, '\nalpha_r=%.6f not achieved hence look up all previoulsy visited solutions.\n', alphar);
                                              
        [backtrackcalls, BTx1, BTx1phat, BTdelta, BTflag] = backtrack(orchandle, ...
            problemparam, solverparam, mk, lastRAipseed, ...
            alphar, logfid);        
        if BTflag==1  % An alpha_r-feasible solution was found in the traceback. 
            x1=BTx1;
            x1phat=BTx1phat;
            delta=BTdelta;
        end
        ncalls=ncalls + backtrackcalls;
        fprintf(logfid, '\n===== BACKTRACK ENDED=====\n\n');                                           
    end

	% PRINT TO LOG
    fprintf(logfid, '\nncalls so far = %d, xbest.ghat = %0 .12f, xbest = [', ncalls, x1.fn);
    fprintf(logfid, '%d ', x1.x);	
    fprintf(logfid, ']\noutput iseed=[');
    fprintf(logfid, '%d ', iseed);	
    fprintf(logfid, ']\n');
    fprintf(logfid, 'secMeas = [');
    fprintf(logfid, '%.12f ', x1.constraint);
    fprintf(logfid, ']\nsecMeasCov = [');
    fprintf(logfid, '%.12f ', x1.ConstraintCov);
    fprintf(logfid, ']\n');
        
    % Print all solutions in stack as column vectors
    for j=1:size(trajectory,1) 
        for i=1:size(trajectory,2); % solutions are stored as columns in stack
            fprintf(logfid, '%d\t', trajectory(j, i));
        end
        fprintf(logfid,'\n');
    end

    %%=====================================================================
    %% THIS PART OF THE CODE IS APPLICABLE ONLY FOR THE FLOWLINE PROBLEM! 
    %% IT WAS INTRODUCED TO PRINT TRUE FUNCTION VALUES TO THE LOG FILES. 
if  TSFtrueFlag == 1 
        htruestar=TSFtrue(problemparam, [6;7;7;12]);
        gtrueplus  = zeros(4,1);
        gtrueminus = zeros(4,1);
        htrueplus  = zeros(4,1);
        htrueminus = zeros(4,1);
        gtrue      = x1.fn;  %deterministic objective
        htrue      = -TSFtrue(problemparam, x1.x) + htruestar;
        ghtrueval  = [gtrue;htrue];
        if htrue<=0
         F='F';
         localmin=1; % feasible
         for p=1:4
            x1.x(p)=x1.x(p)-1;
            gtrueminus(p)=x1.x(1)+x1.x(2)+x1.x(3);
            if x1.x(4)>0
                htrueminus(p)=-TSFtrue(problemparam, x1.x)+htruestar;
            else
                htrueminus(p)=99999;
            end

             x1.x(p)=x1.x(p)+2;
             gtrueplus(p)=x1.x(1)+x1.x(2)+x1.x(3);
             if x1.x(4)<20
               %  htrueplus(p)=-TSFtrue(problemparam, x1.x)+htruestar;                        
             else
                % htrueplus(p)=99999;
             end

             x1.x(p)=x1.x(p)-1;

             if (gtrueplus(p)<gtrue && htrueplus(p)<=0) || (gtrueminus(p)<gtrue && htrueminus(p)<=0)
                 localmin=0;
                 break
             end
         end
         if localmin==1
             if gtrue==20
                 Optm='G';
             else
                 Optm='L';
            end
           FO=Optm;
        else
             Optm='No';
             FO=F;
         end
     else
         F='I';
         Optm='NA';
         FO=F;
     end

   %%======================================================================

     PRINT TO LOG: Compate previous fn and current fn
     fprintf(logfid, 'x1fn=%.6f\n', x1.fn);

     PRINT TO LOG 2: bestsols.txt
     fprintf(logfid2, '% 3d% 10.2f% 7d% 10d [', k, delta, mk, ncalls);
     fprintf(logfid2, '% 8d', x1.x);	
     fprintf(logfid2, ' ]% 14.6f% 9.4f% 15.9f% 17.4f|% 14.6f% 14.9f% 9s% 7s\n', x1.fn, sqrt(x1.FnVar), x1.constraint(1), sqrt(x1.ConstraintCov(1)), ghtrueval(1), ghtrueval(2), F, Optm);
 end  


    %% AD-HOC EARLY TERMINITATION CRITERION:
    %% Exit R-SPLINE if same solution obtained three times in a row. 
    %% NOT TO BE USED. 
%    if isequal(x1,xk)==1
%	countsame=countsame+1;
%    else
%	countsame=0;
%    end
%    if countsame >= 2
% 	break;
%    end

    % Exit R-SPLINE if if budget is exceeded
    if k>=kmax || ncalls>=budget  
        FO ='F';        %BG modification 
        break;
    end
    k=k+1;
end

fprintf(logfid, '\n\n========= R-SPLINE ended after %d iterations=========\n\n', k);
fprintf(logfid, 'Final random-number seed = [');
fprintf(logfid, '%d ',iseed);	
fprintf(logfid, ']\n');
fprintf(logfid, 'Final desision-variable = [');
fprintf(logfid, '%d ', x1.x);	
fprintf(logfid, ']\n');
fprintf(logfid, 'Final objective function value = %0.12f\n', x1.fn);
fprintf(logfid, 'Final constraint function value = [');
fprintf(logfid, '%.12f ', x1.constraint);
fprintf(logfid, ']\n');
fprintf(logfid, 'R-SPLINE ended after %d calls\n', ncalls);
    
end



%% ========================================================================
%  SPLINE
function [ncalls, stackflag, NE_best, iseed]=SPLINE(orchandle, ...
    problemparam, solverparam, x0, mk, bk, k, iseed, logfid)
% =========================================================================

% Declare global variables
global trajectory;
global AirportSecurityFlag
global mobs
% Problem parameters
% id = int32(problemparam(2));	
% nseeds = int32(problemparam(3));
% nsecMeas = int32(problemparam(4));

% Solver paramater
ftol = solverparam(18);

% Initialize SPLINE parameters
NE_best     =struct('x', {}, 'fn', {}, 'FnVar', {}, 'FnGrad', {}, 'FnGradCov', {}, 'constraint', {}, 'ConstraintCov', {}, 'ConstraintGrad', {}, 'ConstraintGradCov', {});
stackflag   =0;  % Set to 1 if all solutions in stack are deemed infeasible at cuurent sampling level mk, else 0.
ncalls      =0;
iseedk      =iseed;
NE_best(1).x   =x0;

fprintf(logfid, 'Initial solution = [');
fprintf(logfid, '%d ', NE_best(1).x);
fprintf(logfid, ']\n');

% Find a sample-path feasible solution in trajectory  
% Move this section of the code into function RSPLINE
stack_ctr=size(trajectory,2); %Each column is a new vector
while (1)
    iseed=iseedk;
    [flag, NE_best, iseed] = orchandle(problemparam, NE_best.x, mk, iseed);
    if flag==0
        if AirportSecurityFlag == 1
        ncalls = ncalls + mobs; %mobs: Number of observations generated
        else
        ncalls = ncalls + mk; % Count mk only if x\in\mathbb{X}.
        end
                
        epsilon=min(max(sqrt(NE_best.ConstraintCov*mk)/mk^0.45, -NE_best.constraint), sqrt(NE_best.ConstraintCov*mk)/mk^0.1);
        delta=(log(sqrt(NE_best.ConstraintCov*mk))-log(epsilon))/log(mk); %delta is a vector of size nsecMeas
        delta(isnan(delta))=solverparam(5);
        if sum(NE_best.constraint > sqrt(mk*NE_best.ConstraintCov)./mk.^delta)>0
            flag=1;
        end
    end
 
    if flag==0 
       if k==1
           
           
           PUSH(NE_best.x);
       end
       break
    end
    % Else initial solution is infeasible
    fprintf(logfid, 'Infeasible initial solution\n');
    
    stack_ctr=stack_ctr-1;
    if stack_ctr<=0  % Stack is empty
        stackflag=1; % Denotes that stack is empty! All visited solutions are infeasible at current sampling level mk
        return       % Returns control to R-SPLINE, which in turn returns control to cgRSPLINE 
    end              % so that local search may be restarted.  

    %trajectory(:,stackctr)=[]; %this would delete x at location stackctr
    NE_best.x=trajectory(:,stack_ctr); %tra::之前找過的可行解
    fprintf(logfid, '[');
    fprintf(logfid, '%d ', NE_best.x);
    fprintf(logfid, ']: ');
end

% Save feasible initial solution to xinit
xinit=NE_best;

fprintf(logfid, '\t===== BEGIN SPLINE LOOP ===\n');
for i=1:bk
	fprintf(logfid, '\t\t=== bk = %d ===\n', i);
		
    iseed=iseedk;   
	[SPLIncalls, SPLI_best, iseed] = SPLI(orchandle, problemparam, ...
        solverparam, NE_best, mk, iseed, logfid);
	
		fprintf(logfid, '\t\t\t==== SPLI ends\n\t\t\tSPLI.x = [');	
	    fprintf(logfid, '%d ', SPLI_best.x);
		fprintf(logfid, ']\n\t\t\tSPLI.fn = %.12f\n', SPLI_best.fn);
		fprintf(logfid, '\t\t\topseed = [');	
	    fprintf(logfid, '%d ', iseed);
		fprintf(logfid, ']\n');
        fprintf(logfid, '\t\t\tsecMeas = [');
	    fprintf(logfid, '%.12f ', double(SPLI_best.constraint));
		fprintf(logfid, ']\n');
	disp(['xold:' num2str(SPLI_best.x(1)) ',' num2str(SPLI_best.x(2)) ',SPLIncalls=' num2str(SPLIncalls)  ])	
    
    iseed=iseedk;
	[NEncalls, NE_best, iseed] = NE(orchandle, problemparam, ...
        solverparam, SPLI_best, mk, iseed, logfid);

        fprintf(logfid, '\t\t\t==== NE ends\n\t\t\tNE_best.x = [');	
	    fprintf(logfid, '%d ', NE_best.x);
		fprintf(logfid, ']\n\t\t\tNE.fn = %.12f\n', NE_best.fn);
		fprintf(logfid, '\t\t\topseed = [');	
	    fprintf(logfid, '%d ', iseed);
		fprintf(logfid, ']\n');
        fprintf(logfid, '\t\t\tsecMeas = [');
	    fprintf(logfid, '%.12f ', double(NE_best.constraint));
		fprintf(logfid, ']\n');
   disp(['xnew:' num2str(NE_best.x(1)) ',' num2str(NE_best.x(2)) ',NEncalls=' num2str(NEncalls)  ])

    ncalls = ncalls + SPLIncalls + NEncalls;
    fprintf(logfid, '\n\n\t\tSPLINE ncalls = prev ncalls + SPLI calls + NE calls = %d\n', ncalls);
        
    if  SPLI_best.x == NE_best.x 
    	fprintf(logfid, '\n\t\tSPLINE ended at bk=%d since NE and SPLI returned the same solution\n\n', i);
    	break
    end
end    
	
% Starting solution is worse than the new solution by at most amount ftol
if xinit.fn <= NE_best.fn + ftol 
	NE_best=xinit;
end

end


%% ========================================================================
%  NE
function [ncalls, NE_best, iseed] = NE(orchandle, problemparam, ...
    solverparam, SPLI_best, mk, iseed, logfid)
global AirportSecurityFlag
global mobs
%  ========================================================================

% Problem parameters
id = int32(problemparam(2));
% nseeds = int32(problemparam(3));
% nsecMeas = int32(problemparam(4));

% Solver parameters
ftol = solverparam(18);

% Initialize NE variables
ncalls    = 0;
iseedk    = iseed;
NE_best   = SPLI_best;
y2        = NE_best.fn;
ixquad    = struct('x', {}, 'fn', {}, 'FnVar', {}, 'FnGrad', {}, 'FnGradCov', {}, 'constraint', {}, 'ConstraintCov', {}, 'ConstraintGrad', {}, 'ConstraintGradCov', {}); 
ixquad(1).x = zeros(id, 1);


fprintf(logfid, '\n\t\t\t===== NE BEGINS =====\n');
fprintf(logfid, '\t\t\tfn at center = %.12f, center = [', y2);
fprintf(logfid, '%d ', NE_best.x);
fprintf(logfid, ']\n');

for i=1:id
	count=1;

    SPLI_best(1).x(i)=SPLI_best(1).x(i)+1;
    iseed=iseedk;
    [flag, SPLI_best, iseed] = orchandle(problemparam, SPLI_best.x, mk, iseed); %#ok<ASGLU>
    if flag==0
        if AirportSecurityFlag == 1
            ncalls = ncalls + mobs; %mobs: Number of observations generated
        else
            ncalls = ncalls + mk;
        end
        epsilon=min(max(sqrt(SPLI_best.ConstraintCov*mk)/mk^0.45, -SPLI_best.constraint), sqrt(SPLI_best.ConstraintCov*mk)/mk^0.1);
        delta=(log(sqrt(SPLI_best.ConstraintCov*mk))-log(epsilon))/log(mk); %delta is a vector of size nsecMeas
        delta(isnan(delta))=solverparam(5);
        if sum(SPLI_best.constraint > sqrt(mk*SPLI_best.ConstraintCov)./mk.^delta)>0
            flag=1;
        end
    end
    if flag == 0  %if feasible, update ncalls and call fn as y2
		y1=SPLI_best.fn; 
		count=count+1;             %count keeps a track of number of feasible points along a direction
        if SPLI_best.fn <= NE_best.fn + ftol % if superior, update current best solution
			NE_best=SPLI_best;
            PUSH(NE_best.x);
        end
    end

    SPLI_best(1).x(i)=SPLI_best(1).x(i)-2;
    iseed=iseedk;
    [flag, SPLI_best, iseed] = orchandle(problemparam, SPLI_best.x, mk, iseed); %#ok<ASGLU>
    if flag==0
        if AirportSecurityFlag == 1
            ncalls = ncalls + mobs; %mobs: Number of observations generated
        else
            ncalls = ncalls + mk;
        end
        epsilon=min(max(sqrt(SPLI_best.ConstraintCov*mk)/mk^0.45, -SPLI_best.constraint), sqrt(SPLI_best.ConstraintCov*mk)/mk^0.1);
        delta=(log(sqrt(SPLI_best.ConstraintCov*mk))-log(epsilon))/log(mk); %delta is a vector of size nsecMeas
        delta(isnan(delta))=solverparam(5);
        if sum(SPLI_best.constraint > sqrt(mk*SPLI_best.ConstraintCov)./mk.^delta)>0
            flag=1;
        end
    end
    if flag == 0
		y3=SPLI_best.fn; 
		count=count+1;
        if SPLI_best.fn <= NE_best.fn + ftol
			NE_best=SPLI_best;
            PUSH(NE_best.x);
        end 
    end

	SPLI_best(1).x(i)=SPLI_best(1).x(i)+1;
	xqnew=SPLI_best.x(i);

    % Quadratic search
    if count==3 
		a = (y1+y3)/2.0 - y2;
		b = (y1-y3)/2.0;
        if a-ftol > 0
			xqnew = int32(SPLI_best.x(i) - (b / (a + a)));
        end
		fprintf(logfid, '\t\t\t\ti = %d, a = %.12f, b = %.12f, xqnew = %.12f\n', i, a, b, xqnew);
		fprintf(logfid, '\t\t\t\ty2 = %.12f, y1 = %.12f, y3 = %.12f\n', y2, y1, y3);
    end
    if  abs(xqnew) < 2147483646.0 %2^31-2
		ixquad(1).x(i) = xqnew;
    end
	fprintf(logfid, '\t\t\t\txold[%d] = %d, ixquad[%d] = %d\n', i, SPLI_best.x(i), i, ixquad.x(i));
end
		
%Call oracle at ixquad. Update NE_best.
iseed=iseedk;
[flag, ixquad, iseed] = orchandle(problemparam, ixquad.x, mk, iseed);
if flag==0
      if AirportSecurityFlag == 1
          ncalls = ncalls + mobs; %mobs: Number of observations generated
      else
          ncalls = ncalls + mk;
      end
    epsilon=min(max(sqrt(ixquad.ConstraintCov*mk)/mk^0.45, -ixquad.constraint), sqrt(ixquad.ConstraintCov*mk)/mk^0.1);
    delta=(log(sqrt(ixquad.ConstraintCov*mk))-log(epsilon))/log(mk); %delta is a vector of size nsecMeas
    delta(isnan(delta))=solverparam(5);
    if sum(ixquad.constraint > sqrt(mk*ixquad.ConstraintCov)./mk.^delta)>0
         flag=1;
    end
end
if flag==0
    if ixquad.fn <= NE_best.fn + ftol
		NE_best=ixquad;
        PUSH(NE_best.x);
    end
end

fprintf(logfid, '\t\t\tixquad.fn = %.12f, ixquad = [', ixquad.fn);
fprintf(logfid, '%d ', ixquad.x);
fprintf(logfid, ']\n');
fprintf(logfid, '\t\t\tixquad.secMeas = [');
fprintf(logfid, '%.12f ', double(ixquad.constraint));
fprintf(logfid, ']\n');
fprintf(logfid, '\t\t\topseed = [');
fprintf(logfid, '%d ', iseed);
fprintf(logfid, ']\n');

end


%% ========================================================================
%  SPLI
function [ncalls, SPLI_best, iseed] = SPLI(orchandle, problemparam, ...
    solverparam, SPLI_best, mk, iseed, logfid)
global AirportSecurityFlag
global mobs
%  ========================================================================

% Problem parameters
id = int32(problemparam(2));
% nseeds = int32(param(3));
% nsecMeas = int32(param(4));

%Solver parameter
ftol = solverparam(18);


% Initialize SPLI variables
ix1 =struct('x', {}, 'fn', {}, 'FnVar', {}, 'FnGrad', {}, 'FnGradCov', {}, 'constraint', {}, 'ConstraintCov', {}, 'ConstraintGrad', {}, 'ConstraintGradCov', {});
imax=100;
jmax=5;
ncalls=0;
s0=1.1;
c=1.1;
iseedk=iseed;

fprintf(logfid, '\t\t\tBEGIN SPLI ===\n');
for j=0:jmax

    fprintf(logfid, '\t\t\t\t=== j = %d ===\n', j);
    fprintf(logfid, '\t\t\t\ti/p seed to PERTURB = [');
    fprintf(logfid, '%d ', iseed);
    fprintf(logfid, ']');
    
    x1=PERTURB(SPLI_best.x, iseed, logfid);
			
    fprintf(logfid, '\t\t\t\tPertrubed value of xbest.x = [');
    fprintf(logfid, '%0.12f ', x1);
    fprintf(logfid, ']\n');

	fprintf(logfid, '\t\t\t\tPLI begins ===\n');		

    iseed=iseedk;
    [PLIncalls, fbar, gamma, npoints, pli_best, iseed] = PLI(orchandle, ...
        problemparam, solverparam, x1, SPLI_best.x, mk, iseed, logfid);
   
    fprintf(logfid, '\t\t\t\t=== PLI ends\n');
    fprintf(logfid, '\t\t\t\tgbar = %.12f\n\t\t\t\tgamma = [',fbar);
    fprintf(logfid, '%.12f ', gamma);
    fprintf(logfid, ']\n\t\t\t\txbest.fn = %.12f\n\t\t\t\txbest = [',pli_best.fn);
    fprintf(logfid, '%d ', pli_best.x);
    fprintf(logfid, ']\n\t\t\t\tsecMeasbest=[');
    fprintf(logfid, '%.12f ', double(pli_best.constraint));
    fprintf(logfid, ']\n');
    fprintf(logfid,'\t\t\t\tnpoints=%d\n', npoints);

	ncalls = ncalls + PLIncalls;	

	% Regardless of whether npoints=id+1 or not, update current best
    if  pli_best.fn + ftol < SPLI_best.fn && npoints>0
        SPLI_best=pli_best;
        PUSH(SPLI_best.x);
    end
		
    if npoints < id+1 % Return control to SPLINE if PLI could not identify 
                      % id+1 feasible points 
        return 
    end
	
    % Perform line search only if PLI finds id+1 feasible points
	glength=norm(gamma);
	fprintf(logfid, '\t\t\t\tglength = %.12f\n', glength);		
    if glength + ftol <= 0
		return
    end
    
    x0=SPLI_best.x;
    gamma=gamma/glength;
    
    for i=0:imax
        s = s0 * c^i;
        ix1(1).x=zeros(id,1);
        for k=1:id
			ix1(1).x(k)=floor(x0(k)-s*gamma(k)+0.5);
        end
        iseed=iseedk;
        [flag, ix1, iseed] = orchandle(problemparam, ix1.x, mk, iseed);
        if flag==0
           if AirportSecurityFlag == 1
            ncalls = ncalls + mobs; %mobs: Number of observations generated
           else
            ncalls = ncalls + mk;
           end
          
            epsilon=min(max(sqrt(ix1.ConstraintCov*mk)/mk^0.45, -ix1.constraint), sqrt(ix1.ConstraintCov*mk)/mk^0.1);
            delta=(log(sqrt(ix1.ConstraintCov*mk))-log(epsilon))/log(mk); %delta is a vector of size nsecMeas
            delta(isnan(delta))=solverparam(5);
            if sum(ix1.constraint > sqrt(mk*ix1.ConstraintCov)./mk.^delta)>0
                flag=1;
            end
        end
        
        %PRINT
        fprintf(logfid, '\t\t\t\t\t== i = %d ==\n\t\t\t\t\ts = %.12f, ix1 = [', i, s);
        fprintf(logfid, '%d ', ix1.x);
        fprintf(logfid, ']\n\t\t\t\t\tix1.fn = %.12f, flag2 = %d\n', ix1.fn, flag);
        fprintf(logfid, '\t\t\t\t\tsecMeas=[');
        fprintf(logfid, '%.12f ', double(ix1.constraint));
        fprintf(logfid, ']\n');
        %END PRINT
	
        if flag~=0 % Return conntrol to SPLINE if line search encounters
                    % an infeasible pt
            return	% ix1 is infeasible
        end
        
        if ix1.fn >= SPLI_best.fn + ftol && i <= 2
            return  % Return control to SPLINE if a worse solution is encountered
                    % right away
        end
        if ix1.fn >= SPLI_best.fn + ftol
            break   % If a worse solution is found later, simply perform PLI again. 
        end
        
        % Update xbest after every step in the line search
		SPLI_best=ix1;
        PUSH(SPLI_best.x);
    end    
end
end



%% ========================================================================
%  PLI
function [ncalls, fbar, gamma, npoints, PLI_best, iseed] = PLI(orchandle, ...
    problemparam, solverparam, x, xbest, mk, iseed, logfid)
global AirportSecurityFlag
global mobs
% Note: PLI does not estimate prob. of feasibility.
% =========================================================================

% Problem parameters
id = int32(problemparam(2));
%nseeds=int32(problemparam(3));
nsecMeas = int32(problemparam(4));

%Solver parameter
ftol = solverparam(18);

% Initialize variables for PLI
PLI_best    =struct('x', {}, 'fn', {}, 'FnVar', {}, 'FnGrad', {}, 'FnGradCov', {}, 'constraint', {}, 'ConstraintCov', {}, 'ConstraintGrad', {}, 'ConstraintGradCov', {});
x0          =struct('x', {}, 'fn', {}, 'FnVar', {}, 'FnGrad', {}, 'FnGradCov', {}, 'constraint', {}, 'ConstraintCov', {}, 'ConstraintGrad', {}, 'ConstraintGradCov', {});
gamma       =zeros(1,id);
npoints     =0;
ncalls      =0;
iseedk      =iseed;
PLI_best(1).x = xbest;
x0(1).x     =floor(x);
strange     =3.145962987654;

z=x-x0.x;  
z=[1;z;0]; % z is a column vector 
[~, p]=sort(z, 'descend');
w=zeros(id+1,1); % w is a column vector
for i=1:id+1		
	w(i)=z(p(i))-z(p(i+1));
end	
wsum=0;
fbar=0;

fprintf(logfid, '\t\t\t\t\tp z w\n');
for i=1:id+1
    fprintf(logfid, '\t\t\t\t\t%d %0.15f %0.15f\n', p(i), z(i), w(i));
end
fprintf(logfid, '\n\t\t\t\t\ti=1, ipseed = [');
fprintf(logfid, '%d ', iseed);
fprintf(logfid, '], mk=%d\n', mk);	
fprintf(logfid, '\t\t\t\t\tx0 = [');
fprintf(logfid, '%d ', x0.x);
fprintf(logfid, ']\n');

% Call oracle at x0
[flag, x0, iseed]=orchandle(problemparam, x0.x, mk, iseed);	
if flag == 0
    if AirportSecurityFlag == 1
       ncalls = ncalls + mobs; %mobs: Number of observations generated
    else
       ncalls = ncalls + mk;
    end 
    epsilon=min(max(sqrt(x0.ConstraintCov*mk)/mk^0.45, -x0.constraint), sqrt(x0.ConstraintCov*mk)/mk^0.1);
    delta=(log(sqrt(x0.ConstraintCov*mk))-log(epsilon))/log(mk); %delta is a vector of size nsecMeas
    delta(isnan(delta))=solverparam(5);
    if sum(x0.constraint > sqrt(mk*x0.ConstraintCov)./mk.^delta)>0
        flag=1;
    end
end
if flag == 0
    npoints=npoints+1;
    wsum = wsum + w(1);
    fbar = fbar + w(1)*x0.fn;
    ghatold = x0.fn;
    PLI_best=x0;

    %PRINT TO LOG
    fprintf(logfid, '\t\t\t\t\t secMeas = [');
    for i=1:nsecMeas
        fprintf(logfid, '%0.12f', double(x0.constraint(i)));
    end
    fprintf(logfid, ']\n');
    fprintf(logfid, '\t\t\t\t\tx0.fn = %.12f, w = %.12f, ghatold = %.12f\n\t\t\t\t\topseed = [', x0.fn, w(1), ghatold);
    fprintf(logfid, '%d ', iseed);
    fprintf(logfid, '], gbar=%.6f\n', fbar);	
    %END PRINT			
else
    ghatold = 0;
    PLI_best(1).fn = strange;
end

% Call oracle at the other id points that form the simplex
for i=2:id+1  

    x0.x(p(i)-1)=x0.x(p(i)-1)+1;
    iseed=iseedk;
    fprintf(logfid, '\n\t\t\t\t\ti = %d, ipseed = [', i);
    fprintf(logfid, '%d ', iseed);
    fprintf(logfid, '], mk = %d\n', mk);		
	fprintf(logfid, '\t\t\t\t\tx0 = [');
    fprintf(logfid, '%d ', x0.x);
    fprintf(logfid, ']\n');

	[flag, x0, iseed]=orchandle(problemparam, x0.x, mk, iseed);
    if flag == 0
        if AirportSecurityFlag == 1
            ncalls = ncalls + mobs; %mobs: Number of observations generated
        else
            ncalls = ncalls + mk;
        end
    
        epsilon=min(max(sqrt(x0.ConstraintCov*mk)/mk^0.45, -x0.constraint), sqrt(x0.ConstraintCov*mk)/mk^0.1);
        delta=(log(sqrt(x0.ConstraintCov*mk))-log(epsilon))/log(mk); %delta is a vector of size nsecMeas
        delta(isnan(delta))=solverparam(5);
        if sum(x0.constraint > sqrt(mk*x0.ConstraintCov)./mk.^delta)>0
            flag=1;
        end
    end
    if flag == 0
        npoints=npoints+1;
        wsum = wsum + w(i);
        fbar = fbar + w(i)*x0.fn;
        gamma(p(i)-1) = x0.fn - ghatold;
        ghatold = x0.fn;

        %PRINTS
        fprintf(logfid, '\t\t\t\t\tx0.fn = %.12f, w = %.12f, ghatold = %.12f\n\t\t\t\t\topseed = [', x0.fn, w(i), ghatold);
        fprintf(logfid, '%d ', iseed);
        fprintf(logfid, '], gbar=%0.6f, npoints=%d\n', fbar, npoints);
        %END PRINTS		

        if PLI_best.fn == strange || x0.fn < PLI_best.fn + ftol
            PLI_best=x0;
        end

        %PRINTS
        fprintf(logfid, '\t\t\t\t\txbest.fn = %.12f\n\t\t\t\t\txbest = [', PLI_best.fn);
        fprintf(logfid, '%d ', PLI_best.x);
        fprintf(logfid, ']\n');
        fprintf(logfid, '\t\t\t\t\tsecMeas = [');
        fprintf(logfid, '%.12f ', PLI_best.constraint);
        fprintf(logfid, ']\n\n');
        %END PRINTS
    end
end

if wsum > ftol
    fbar = fbar/wsum;	
end

fprintf(logfid, '\n\t\t\t\t\tgbar = %.12f, npoints = %d\n', fbar, npoints);
end



%% ========================================================================
%  PERTURB
function xpert = PERTURB(x, iseed, logfid)

% =========================================================================

id=length(x);
xpert=zeros(id,1);  % Column vector
seed=iseed(1);

fprintf(logfid, '\t\t\t\tPERTURB begins ==== \n');
for i=1:id
    fprintf(logfid, '\t\t\t\t\t ipseed = [');
    fprintf(logfid, '%d ', seed);
    fprintf(logfid, '], ');
        
	[seed, u] = u16807d(seed);
    xpert(i) = x(i) + .3*(u - 0.5);

    fprintf(logfid, 'u = %0.12f, xpert[%d] = %.12f\n', u, i, xpert(i));

end
fprintf(logfid, '\t\t\t\t==== PERTURB ends\n');
end




%% ========================================================================
% 
function [ncalls, BTxnew, BTxnewphat, BTdelta, BTflag] = backtrack(orchandle, ...
    problemparam, solverparam, mk, iseed, alphar, logfid)
global AirportSecurityFlag
global mobs
% Searches trajectory for an alphar-feasible solution.
% If one isn't found, set BTflag to 0. (BTflag=1 if alpha_r feas sol is
% found)
% =========================================================================

% Declare global variables
global trajectory

% Problem parameters
id=int32(problemparam(2));
% nseeds=int32(problemparam(3));
% nsecMeas=int32(problemparam(4));

% Initialize backtrack varaibles
ncalls   =0;
iseedk   =iseed;
stack_ctr =size(trajectory,2); %Each column is a new vector
BTflag   =0;  % Alpha_r feasible solution not found

fprintf(logfid, 'Trajectory:\n');
for i=1:id
    for j=1:stack_ctr
        fprintf(logfid, '%d ', trajectory(i,j));
    end
    fprintf(logfid, '\n');
    
end

while stack_ctr>0
    %trajectory(:,stack_ctr)=[]; %this would delete x at location stackctr
    iseed=iseedk;
    [flag, BTxnew, iseed] = orchandle(problemparam, trajectory(:,stack_ctr), mk, iseed); %#ok<ASGLU>
    
    if flag == 0
        if AirportSecurityFlag == 1
          ncalls = ncalls + mobs; %mobs: Number of observations generated
        else
          ncalls = ncalls + mk;
        end    % Count mk only if x\in\mathbb{X}. 

        epsilon=min(max(sqrt(BTxnew.ConstraintCov*mk)/mk^0.45, -BTxnew.constraint), sqrt(BTxnew.ConstraintCov*mk)/mk^0.1);
        BTdelta=(log(sqrt(BTxnew.ConstraintCov*mk))-log(epsilon))/log(mk); %delta is a vector of size nsecMeas
        BTdelta(isnan(BTdelta))=solverparam(5);
        if sum(BTxnew.constraint > sqrt(mk*BTxnew.ConstraintCov)./mk.^BTdelta)>0
            flag=1;
        end
        BTxnewphat=prod(tcdf(mk.^(0.5-BTdelta)-BTxnew.constraint./sqrt(BTxnew.ConstraintCov), mk-1));
    
        if flag == 0 && BTxnewphat>=alphar
           % Initial solution is feasible and alphar-feas.
           BTflag=1;
           break
        end
    end
    stack_ctr=stack_ctr-1;
    
end

% PRINT TO LOG
if stack_ctr>0 
    fprintf(logfid, '\n\tBacktrack found a level alpha_r solution, xnew = [');
    fprintf(logfid, '%d ', BTxnew.x);	
    fprintf(logfid, ']\n');
    fprintf(logfid, '\tsecMeas = [');
    fprintf(logfid, '%.12f ', BTxnew.constraint);
    fprintf(logfid, ']\n\tsecMeasCov = [');
    fprintf(logfid, '%.12f ', BTxnew.ConstraintCov);
    fprintf(logfid, ']\n\txbestfn = %d, phat = %.12f\n', BTxnew.fn, BTxnewphat);
else
    fprintf(logfid, '\n\tBacktrack did not yield a level alpha_r solution.\n');
    % In which case xnew being returned is the first solution to enter the
    % stack
end

end





%% ========================================================================
%%PUSH
function PUSH(x)
%  ========================================================================
global trajectory;
x=reshape(x,[numel(x),1]);

% If trajectory is empty, set trajectory=x
if numel(trajectory)==0
    trajectory=x;
    return
end

% Search for x in trajectory. Append x to trajectory if not found.
if isempty(find(ismember(trajectory',x','rows'),1))        %BG modification 
    trajectory=[trajectory x];
end
disp('trajectory')
disp( trajectory)
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
