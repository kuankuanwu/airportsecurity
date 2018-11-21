function [flag, xstruct, iseed] = ThreeStageFlowline(problemparam, x, m, iseed)
	
% =========================================================================
%     input:
%       param = 
%         param(1) = problem ID, not used
%         param(2) = problem dimension = 4
%         param(3) = nseeds = 3
%         param(4) = nSecMeas = 1 
%         param(5) = warm-up time
%         param(6) = simulation end time
%         param(7) = total service rate
%         param(8) = total buffer space available
%       ix    =
%         integer solution, ( s1, s2, s3, b2 )
%       m     = sample size
%       iseed = 
%         iseed(1), iseed(2), iseed(3) = one seed per server
%     output:
%       iflag1 = 0 implies that the model parameters are feasible
%              = 1 implies that the model parameters are infeasible
%       iflag2 = 0 implies that ix is feasible
%              = 1 implies that ix is infeasible
%       fn     =      ybar (defined only if flag1 = 0 and flag2 = 0)	
%       constraint = matrix (size = 1 X nsecMeas) of estimates of 
%                   constraint functions 
%
% =========================================================================
% Example: 
% param = [123 4 3 1 50 1000 20 20]
% x=[6;7;7;12]
% m=100
% iseed=[23942;209482;348]
% 
% =========================================================================


% Problem parameters
% id     = int32(problemparam(2));
% nseeds = int32(problemparam(3));
% nsecMeas = int32(problemparam(4));
warmuptime = problemparam(5);
totaltime = problemparam(6);
ratetotal = problemparam(7);
buffertotal = problemparam(8);

%Initialize variables
xstruct=struct('x', {}, 'fn', {}, 'FnVar', {}, 'FnGrad', {}, 'FnGradCov', {}, 'constraint', {}, 'ConstraintCov', {}, 'ConstraintGrad', {}, 'ConstraintGradCov', {});
xstruct(1).x=x;

% Check for feasibility
flag=0;
if warmuptime < 0 || totaltime <= warmuptime || buffertotal < 0
    flag=1;
    %fprintf(1, 'INVALID MODEL PARAMETERS!!\n');
	return
end
rate1 = x(1);
rate2 = x(2);
rate3 = x(3); % ratetotal - rate1 - rate2 + 0.5;
buf2 =  x(4);  % x(3);
buf3 = buffertotal - buf2 ;
if sum(x<=0)>0 || sum(x(1:3)>ratetotal)>0 || x(4)>buffertotal-1 
    flag = 1;
    %fprintf(1, 'INFEASIBLE x!\n');
    return
end

sum1 = 0;
sum2 = 0;
for i=1:m	%m days of simulation of the system
    icount = 0;     
	%fprintf(1, 'rep = %d\n', i);
		
	%simulate: return throughput / (total time - warmup time)
	%initialize clock, state, and event calendar

	timebig = totaltime + 1;
	thruput = 0;
	%clock = 0;
	nbuf2 = 0;
	nbuf3 = 0;
		
    [iseed(1), u1] = u16807d(iseed(1));
	tend1 = -log(1 - u1)/rate1;
	tend2 = timebig;
	tend3 = timebig;
        %fprintf(1, 'seed1, seed2, seed3 =%d, %d, %d\n', iseed(1), iseed(2), iseed(3));
        %fprintf(1, 'tend1, tend2, tend3 = %.6f, %.6f, %.6f\n', tend1, tend2, tend3);

    while (1)		
		%//get next-event time and number
		timemin = totaltime;
		nevent = 4;
        if tend1 <= timemin
			timemin = tend1;
			nevent = 1;
        end
        if tend2 <= timemin
			timemin = tend2;
			nevent =2;
        end
        if tend3 <= timemin
            timemin = tend3;
			nevent =3;
        end
		clock = timemin;
			
		%execute the next event
        if nevent == 1
            %//server 1: end of service
            if nbuf2 == buf2 
                tend1 = timebig; 	%//buffer of server 2 is full
            else
				nbuf2 = nbuf2 + 1;
				[iseed(1), u1] = u16807d(iseed(1));
				tend1 = clock -log(1 - u1)/rate1;
                if nbuf2 == 1
					[iseed(2), u2] = u16807d(iseed(2));
					tend2 = clock -log(1 - u2)/rate2;
                end
            end
        elseif nevent == 2
			%//server 2: end of service
            if nbuf3 == buf3
                tend2 = timebig;
            else
				nbuf2 = nbuf2 - 1;
				nbuf3 = nbuf3 + 1;
                if tend1 == timebig
					nbuf2 = nbuf2 + 1;
    				[iseed(1), u1] = u16807d(iseed(1));
    				tend1 = clock -log(1 - u1)/rate1;
                end
                if nbuf2 > 0
                    [iseed(2), u2] = u16807d(iseed(2));
                    tend2 = clock -log(1 - u2)/rate2;
                else
                    tend2 = timebig;
                end
                if nbuf3 == 1
                    [iseed(3), u3] = u16807d(iseed(3));
                    tend3 = clock -log(1 - u3)/rate3;
                end
            end
        elseif nevent == 3
			%//server 3: end of service
            if clock >= warmuptime
                thruput = thruput + 1;
            end
			nbuf3 = nbuf3 - 1;
            if nbuf2 > 0 && tend2 == timebig
				nbuf2 = nbuf2 - 1;
				nbuf3 = nbuf3 + 1;
                if nbuf2 > 0
                    [iseed(2), u2] = u16807d(iseed(2));
                    tend2 = clock -log(1 - u2)/rate2;
                end
            end
            if nbuf2 < buf2 && tend1 == timebig
				nbuf2 = nbuf2 + 1;
				[iseed(1), u1] = u16807d(iseed(1));
				tend1 = clock -log(1 - u1)/rate1;
                if nbuf2 == 1
                    [iseed(2), u2] = u16807d(iseed(2));
                    tend2 = clock -log(1 - u2)/rate2;
                end
            end
			tend3 = timebig;
            if nbuf3 > 0
				[iseed(3), u3] = u16807d(iseed(3));
				tend3 = clock -log(1 - u3)/rate3;
            end
        elseif nevent == 4
            %//event: end of simulation
			%// minimizing: return the negative of throughput
			sum1 = sum1 + (-thruput / (totaltime - warmuptime));
			sum2 = sum2 + (-thruput / (totaltime - warmuptime))^2;
			break
        end
    	icount=icount+1;
            %fprintf(1,'\neven %d finished\n', icount);
            %fprintf(1,'next event = %d, clock = %.12f\n', nevent, clock);
            %fprintf(1,'nbuf2, nbuf3 = %d, %d\n', nbuf2, nbuf3);
            %fprintf(1,'tend1, tend2, tend3 = %0.12f, %.12f, %.12f\n', tend1, tend2, tend3);
            %fprintf(1,'seeds = %d, %d, %d\n", iseed(1), iseed(2), iseed(3));
            %fprintf(1,'throughput = %.12f\n\n", thruput);
    end
end
sum1=sum1/m;
xstruct(1).fn=rate1+rate2+rate3;
xstruct(1).FnVar=0;
xstruct(1).constraint = sum1 + TSFtrue(problemparam, [6;7;7;12]); % Add flowtrue([6 7 7 12]) on LHS to make RHS=0
xstruct(1).ConstraintCov=((sum2/m) - sum1^2)/m;

    %fprintf(1, 'ix = [%d %d %d %d %d]: \n', rate1, rate2, rate3, buf2, buf3);
    %fprintf(1, 'ghat = %.12f, hhat = %0.12f\n', fn, constraint);
    %fprintf(1, 'varhat_hhat=%.6f\n', ConstraintCov);
end



function [iseed,u16807d]=u16807d(iseed)
%..........................................................................
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
u16807d=0;
while (u16807d<=0 || u16807d>=1)
    iseed = mod (iseed * 16807,2147483647);
    u16807d = iseed / 2147483648;
end
end