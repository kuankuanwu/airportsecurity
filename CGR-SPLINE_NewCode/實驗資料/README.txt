===========================================================================
Follow these steps to generate one sample-path solution trajectory of
cgR-SPLINE on the three-stage flowline problem. This is achieved by running
top.m with an input seed vector to the oracle (problemseed) and an input 
seed vector to cgR-SPINE (solverseed). Change both seed values and re-run 
top.m to generate a new sample-path. 
===========================================================================

(1.) For efficiency we suggest generating a mex file for the oracle by 
running this command at the command line:
>> mex TSF.c
If mex ran successfully, MATLAB should return a message that looks 
something like this:

Building with 'Xcode with Clang'.
MEX completed successfully.

See more information, see the MATLAB reference page for building mex files

If you are unable to build a mex file, make the following changes to top.m.  
This should run the cgR-SPLINE with the oracle coded in a .m file. 

problemname      = 'ThreeStageFlowline';

(2.) Run top.m with relevant input parameters. For details see top.m
>> top([1;2;3],[2;3;4;5])


(3.) Running top.m should generate two files
a. TSF.txt          : Logs almost every major algorithmic step in 
                      cgR-SPLINE. 
b. TSF_bestsols.txt : Logs all the solutions returned by SPLINE after 
                      every RA interation in each restart. You will notice 
                      that some restart numbers repeat. This is because the
                      code repeats a restart if it exits a local search
                      before it can locate a sample-path feasible solution.
c. TSF_report.txt   : Contains a table of the current local solution X_r 
                      and the updated incumbent Z_r at the each restart.
d. TSF_final.txt    : Contains a table of just the incumbents, Z_r.
e. TSF_vars.mat     : Generated at the end of the sample-path run once the 
                      total simulation budget is exhausted. 
                      Data from _bestsols, _report, and _final is 
                      stored in TSF_vars.mat. See cgRSPLINE_v2.m for a
                      description of the components of TSF_vars.mat.   

(4.) Files generated from running top.m with problemseed=[1;2;3] and 
     solverseed=[2;3;4;5] are stored in ExampleOutput/