#include<stdio.h>
#include<math.h>
#include "mex.h"
#include "matrix.h"
   
/* Data structure for function u16807d */
typedef struct {
	int iseed;
	long double RV;
} u16807d; 


/* u16807d */ 
u16807d RNG_u16807d(int iseed){

	u16807d next;

	next.RV=0;
	next.iseed=iseed;

	while(next.RV <= 0 || next.RV >= 1){
		next.iseed = (next.iseed * 16807L) % 2147483647;
		next.RV = next.iseed/2147483648.0;
	}
	return next;
}	


/* Generate Poisson random variates using u16807d */
/* CDF-inverse method. */
int poissgen(long double u, int lambda){
	long double s,p;
	int x=0;
	
	p=exp(-lambda);
	s=p;
	while(u>s){
		x=x+1;
		p=p*lambda/x;
		s=s+p;
	}
	return x;
}

void ThreeStageFlowline(double *param, mwSize x1, mwSize x2, mwSize x3, mwSize x4, mwSize m, mwSize iseed1, mwSize iseed2, mwSize iseed3, double *flagseed, double *fn, double *varfn){

/*  THE THREE-STAGE FLOWLINE PROBLEM
    INPUT:
        param = 
        param(0) = problem ID, not used
        param(1) = problem dimension = 4
        param(2) = nseeds = 3
        param(3) = nSecMeas = 1
        param(4) = warm-up time
        param(5) = simulation end time
        param(6) = total service rate
        param(7) = total buffer size

        m     = sample size

    OUTPUT:
        iflag1 = 0 implies that the model parameters are feasible
               = 1 implies that the model parameters are infeasible
        iflag2 = 0 implies that ix is feasible
               = 1 implies that ix is infeasible
        fn     = ybar (defined only if flag1 = 0 and flag2 = 0)
        var    = var(ybar) (defined only if flag1=flag2=0)
*/

        
	int i, icount=0;
	int rate1, rate2, rate3, buf2, buf3;
	int nbuf2, nbuf3, nevent;
	long double warmuptime, totaltime, ratetotal, buffertotal;
	long double sum, sum2, timebig, thruput, clock, tend1, tend2, tend3, timemin;
	int flag;
    u16807d u1, u2, u3;
    
    warmuptime  = param[4];
	totaltime   = param[5];
	ratetotal   = param[6];
	buffertotal = param[7];

	flag=0;
    /* model parameters feasibility check */
	if(warmuptime < 0 || totaltime <= warmuptime || ratetotal < 3  || buffertotal < 0) {flag=1;}

	/* decision-variable feasibility check */
	if(x1<0 || x2<0 || x3<0 || x1>ratetotal || x2>ratetotal || x3>ratetotal || x4>buffertotal-1 ){flag = 1;}
	    
    /* Initialize output */
    flagseed[0]=flag;
    flagseed[1]=iseed1;
    flagseed[2]=iseed2;
    flagseed[3]=iseed3;
    fn[0]      = 0;
    fn[1]      = 0;
    varfn[0]   = 0;
    varfn[1]   = 0;

	/* return if param or x infeasible */
    if (flag == 1) {
		return;
	}

    /* Initialize loop variables */
    rate1 = x1;
	rate2 = x2;
	rate3 = x3; /* ratetotal - rate1 - rate2 + 0.5; */
	buf2  = x4;
	buf3  = buffertotal - buf2;

	u1.iseed = iseed1;
	u2.iseed = iseed2;
	u3.iseed = iseed3;

	sum  = 0;
    sum2 = 0;
	for(i=1;i<=m;i++){	// m days of simulation of the system
	   	icount = 0;     
		
		//printf("iref = %d", i);
		
		//simulate: return throughput / (total time - warmup time)
		//initialize clock, state, and event calendar

		timebig = totaltime + 1;
		thruput = 0;
		clock = 0;
		nbuf2 = 0;
		nbuf3 = 0;
		
		iseed1= u1.iseed;
        u1    = RNG_u16807d(iseed1);
		tend1 = -logl(1.0 - u1.RV)/rate1;
		tend2 = timebig;
		tend3 = timebig;

		//printf("seed1, seed2, seed3 =%d, %d, %d\n", iseed1, iseed2, iseed3);
		//printf("tend1, tend2, tend3 = %Lf, %Lf, %Lf\n", tend1, tend2, tend3);

		while(1){
			
			//get next-event time and number
			timemin = totaltime;
			nevent = 4;
			if(tend1 <= timemin){
				timemin = tend1;
				nevent =1;
			}
			if(tend2 <= timemin){
				timemin = tend2;
				nevent =2;
			}
			if(tend3 <= timemin){
				timemin = tend3;
				nevent =3;
			}
		    clock = timemin;
			
			//execute the next event
			if(nevent == 1){
				//server 1: end of service
				if(nbuf2 == buf2){ tend1 = timebig;}
				else{
					nbuf2 = nbuf2 + 1;
					iseed1=u1.iseed;
			        u1 = RNG_u16807d(iseed1);
					tend1 = clock -logl(1.0 - u1.RV)/rate1;
						if(nbuf2 == 1){
						iseed2=u2.iseed;
						u2 = RNG_u16807d(iseed2);
						tend2 = clock -logl(1.0 - u2.RV)/rate2;
					}
				}
			}
			else if(nevent == 2){
				//server 2: end of service
				if(nbuf3 == buf3) {tend2 = timebig;}
				else{
					nbuf2 = nbuf2 - 1;
					nbuf3 = nbuf3 + 1;
					if(tend1 == timebig){
						nbuf2 = nbuf2 + 1;
						iseed1=u1.iseed;
			        	u1 = RNG_u16807d(iseed1);
						tend1 = clock -logl(1.0 - u1.RV)/rate1;
					}
					if(nbuf2 > 0){
						iseed2=u2.iseed;
			        	u2 = RNG_u16807d(iseed2);
						tend2 = clock -logl(1.0 - u2.RV)/rate2;
					}
					else{ tend2 = timebig;}
					if(nbuf3 == 1){
						iseed3 = u3.iseed;
				        u3 = RNG_u16807d(iseed3);
						tend3 = clock -logl(1.0 - u3.RV)/rate3;
					}
				}
			}
			else if(nevent == 3){
				//server 3: end of service
				if(clock >= warmuptime) {thruput = thruput + 1;}
				nbuf3 = nbuf3 - 1;
				if(nbuf2 > 0 && tend2 == timebig){
					nbuf2 = nbuf2 - 1;
					nbuf3 = nbuf3 + 1;
					if(nbuf2 > 0){
						iseed2=u2.iseed;
						u2 = RNG_u16807d(iseed2);
						tend2 = clock -logl(1.0 - u2.RV)/rate2;
					}
				}
				if(nbuf2 < buf2 && tend1 == timebig){
					nbuf2 = nbuf2 + 1;
					iseed1=u1.iseed;
			        u1 = RNG_u16807d(iseed1);
					tend1 = clock -logl(1.0 - u1.RV)/rate1;
					if(nbuf2 == 1){
						iseed2=u2.iseed;
						u2 = RNG_u16807d(iseed2);
						tend2 = clock -logl(1.0 - u2.RV)/rate2;
					}
				}
				tend3 = timebig;
				if(nbuf3 > 0){
					iseed3=u3.iseed;
					u3 = RNG_u16807d(iseed3);
					tend3 = clock -logl(1.0 - u3.RV)/rate3;
				}
			}			
			else if(nevent == 4){
				// event: end of simulation
				// minimize negative throughput
				sum  = sum  + (-thruput / (totaltime - warmuptime));
                sum2 = sum2 + pow(-thruput / (totaltime - warmuptime),2);
				break;
			}
			icount++;

		    //printf("\neven %d finished\n", icount);
			//printf("next event = %d, clock = %.15Lf\n", nevent, clock);
			//printf("nbuf2, nbuf3 = %d, %d\n", nbuf2, nbuf3);
			//printf("tend1, tend2, tend3 = %0.15Lf, %.15Lf, %.15Lf\n", tend1, tend2, tend3);
			//printf("seeds = %d, %d, %d\n", u1.iseed, u2.iseed, u3.iseed);
			//printf("throughput = %.15Lf\n\n", thruput);

		}
	}
    fn[0]       = x1+x2+x3;
    fn[1]       = sum/m;
    varfn[0]    = 0;
    varfn[1]    = (sum2/m - pow(fn[1],2))/m;
    flagseed[1] = u1.iseed;
    flagseed[2] = u2.iseed;
    flagseed[3] = u3.iseed;

    //printf("orcflow: fn = %Lf, ix = [%d %d %d %d]\n\n", fn[1], x1, x2, x3, x4);


}


/* Gateway Function */
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    double  *param;
    int     x1, x2, x3, x4;
    int     m, inputseed1, inputseed2, inputseed3;            
    double  *flagseed;       
    double  *fn, *varfn;
	
    /* Check for proper number of arguments. */
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("MyToolbox:matrixSum:nrhs","Nine inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:matrixSum:nlhs","Three output required.");
    }

    /* Make a pointer to input data. */
    param = mxGetPr(prhs[0]);
    x1    = mxGetScalar(prhs[1]);
    x2    = mxGetScalar(prhs[2]);
    x3    = mxGetScalar(prhs[3]);
    x4    = mxGetScalar(prhs[4]);
    m     = mxGetScalar(prhs[5]);
    inputseed1 = mxGetScalar(prhs[6]);
    inputseed2 = mxGetScalar(prhs[7]);
    inputseed3 = mxGetScalar(prhs[8]);
         
    /* Make the output matrices. */
    plhs[0] = mxCreateDoubleMatrix(1,4,mxREAL);  /* flagseed */
    plhs[1] = mxCreateDoubleMatrix(2,1,mxREAL);  /* fn       */
    plhs[2] = mxCreateDoubleMatrix(2,1,mxREAL);  /* FnVar    */
    
    /* Make a pointer to the real data in the output matrix. */
    flagseed = mxGetPr(plhs[0]);
    fn       = mxGetPr(plhs[1]);
    varfn    = mxGetPr(plhs[2]);

    /* Call the computational routine. */
    ThreeStageFlowline(param, x1, x2, x3, x4, m, inputseed1, inputseed2, inputseed3, flagseed, fn, varfn);
}