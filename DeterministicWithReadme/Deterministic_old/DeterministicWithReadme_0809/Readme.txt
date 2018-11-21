Readme

Filename description:
Notice: Only main script can be run directly without additonal global variable configuration

A. Main programs
1.	Main_AirportSecurity.m:
	The main script to find the optimal values of (S1,S2,tau) 
	which result in the maximum security level under constraints on both mean waiting time and budget.
	The feasible region of tau is [0,1] and (S1, S2) are nonnegative integers.
	Need input file (Input.txt) with the specified format (Details please refer part C. Input File)
	The output file (Output.txt) recorded the optimization process and the corresponded results.
	
	Subroutines used:
	InitializeAndReadInput.m
	GlobalOptimizaton.m
	CalcR1.m
	CalcR2.m
	
B. Subroutines
1.	InitializeAndReadInput.m
	This function initialize the airport model by reading the input text file.
	Need input file (Input.txt) with the specified format (Details please refer part C. Input File)
2. 	GlobalOptimizaton.m
	Function of global optimization to find the (S1,S2) which can produce the global minimum tau (and maximum security level)
	This function will generate the initial (S1,S2) and then perform the SPLINE optimization.
	It will terminate after it have perform K iteration(K is specified by user in input file mentioned in main scripts)
	The intermediate results will be print on screen and write into output file.

3.	SPLINE.m
	It combines SPLI (Line search) and NE (neigborhood enumeration) to find the (S1,S2) which can produce the local minimum tau (and maximum security level) given the initial (S1,S2)
	For each iteration, SPLI and NE are performed. If the results (S1,S2) from SPLI and NE are the same, it will return this (S1,S2) as the local optimum.

4.	SPLI.m
	Line search to update (S1,S2) in order to produce the lower tau (and higher security level).
	
5.	PLI.m
	Evaluate the gradient using interpolation method for line search procedure.
	
6.	NeighborEnum.m
	Neigborhood enumeration to update (S1,S2) in order to produce the lower tau (and higher security level).
	
7.	GetOptimTau.m
	Note: All optimization method call this function to evaluate the minimum tau by given (S1,S2).
	Finding the minimum tau (and maximum security level) by given (S1,S2).
	First, it will call GetMinTau() to find the minimum tau without the constraint of mean waiting time.
	If the tau can satisfy the constraint of mean waiting time, it will return directly.
	Otherwise, it will call the SolveTau.m to solve the minimum tau with the constraint of mean waiting time.

8.	GetMinTau.m
	Find the minimum tau by given (S1,S2) without the constraint of mean waiting time.

9.	SolveTau.m
	Find the minimum tau by given (S1,S2) with the constraint of mean waiting time.
	
	First, the lower bound of tau is obtained by GetMinTau() and became the input argument of this function. 
	The lower bound of tau will produce mean waiting time > Epslion which is violate the constraint.
	
	To find the upper bound, this function try to find the local minimum of mean waiting time. 
	FindUpperBound_Guess() can produce the reasonable initial guess of local minimum of mean waiting time 
	by finding the tau which can make the traffic rate of Non-selectee lane and Selectee lane the same.
	This tau from FindUpperBound_Guess() usually can satisfy the constraint of mean waiting time.
	If this value cannot satisfy the constraint of mean waiting time, fminsearch() will be applied using this 
	value as the initial value.
	If the local minimum of mean waiting time cannot satisfy the constraint of mean waiting time, it return no solution.
	
	After it obtains the lower bound and the upper bound, it will apply the RegularFalsi_Bisection() to find the 
	minimum tau which satisfy the constraint of mean waiting time.

10.	FindUpperBound_Guess.m
	Find the reasonable initial guess of local minimum of mean waiting time 
	by finding the tau which can make the traffic rate of Non-selectee lane and Selectee lane the same.

11.	RegularFalsi_Bisection.m
	Find the minimum tau which satisfy the constraint of mean waiting time by given the lower bound and the upper bound.
	The lower bound of tau should violate the constraint of mean waiting time.
	The upper bound of tau should satisfy the constraint of mean waiting time.
	
12. GetReminTime.m
	Calculate the difference between the mean waiting time and epsilon.
	Return the negative value if the mean waiting time > epsilon.
	Return the positive value if the mean waiting time < epsilon.

13.	MeanWait.m
	Evaluate the mean waiting time.
	It can call the matrix inverse version (MeanWait_MatrixInverse.m) or chain version (MeanWait_Chain.m).
	The default method is set to matrix inverse to avoid the error propagation.
	
14. MeanWait_MatrixInverse.m
	Evaluate the mean waiting time using matrix inverse method
	
15. MeanWait_Chain.m
	Evaluate the mean waiting time using chain inverse method

16.	GetA.m
	Evaluate the transition rate matrix of non-selectee lane given the current state 
	and the future state of selectee lane.
	
17. CalcR1.m
	Evaluate the risk ratio of non-selectee lane
	
18. CalcR2.m
	Evaluate the risk ratio of selectee lane

C. Input files
	Input.txt: Input file for 
	Line 1: c1 value defined in the thesis
	Line 2: c2 value defined in the thesis
	Line 3: Beta1 value defined in the thesis
	Line 4: Beta2 value defined in the thesis
	Line 5: Budget value defined in the thesis
	Line 6: Epsilon value defined in the thesis
	Line 7: N value defined in the thesis
	Line 8: Lambda value defined in the thesis
	Line 9: mu1 value defined in the thesis
	Line 10: mu2 value defined in the thesis
	Line 11: theta value defined in the thesis
	Line 12: Method for evaluating mean waiting time (1==>MatrixInverse, 2==>Chain)
	Line 13: Number of initial points for optimization
	Line 14: The limit of traversed points (S1,S2) 

