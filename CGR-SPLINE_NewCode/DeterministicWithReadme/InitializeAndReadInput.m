function InitializeAndReadInput(InputFileName)
%  purpose:     To set the parameters of the airport model and set the
%               simulation criteria.
%  input file
%   Input.txt:  Contains the parameters of the systems
%  output file:
%   Output.txt: Log file during optimization. Finally it will record the
%               minimum tau and the maximum security level.
    %(Constraint)
    global c1
    global c2
    global Beta1
    global Beta2
    global Budget
    global Epsilon
    %  c1 = 1;
    %  c2 = 2;
    %  Beta1 = 3;
    %  Beta2 = 15;
    %  Budget = 44;
    %  Epsilon = 10;

    %(Waiting Time)
    global N
    global Lambda
    global mu1
    global mu2
    global theta
    %  N = 100;
    %  Lambda = 41;
    %  mu1 = 8;
    %  mu2 = 6;
    %  theta = 0.0625;
    global MeanWait_Flag
    
    global NInitial
    global seed
    global b
    global d1
    global d2

    global fout
    fid = fopen(InputFileName,'r');
    fout = fopen('Output.txt','w');


    if fid < 0
        disp('Error: File Not Found!');
        return
    end
    disp('Input Constraint Conditions:')
    
    c1 = fscanf(fid,'%f',1);
    disp(['Input c1:',num2str(c1)])
    if c1 <= 0
        disp('Error: c1 should > 0')
        return
    end
    
    c2 = fscanf(fid,'%f',1);
    disp(['Input c2:',num2str(c2)])
    if c2 <= 0
        disp('Error: c2 should > 0')
        return
    end
    if c1 > c2
        disp('Error: c1 should < c2')
        return
    end
    Beta1 = fscanf(fid,'%f',1);
    disp(['Input Beta1:',num2str(Beta1)])
    if Beta1 <= 0
        disp('Error: Beta1 should > 0')
        return
    end
   
    Beta2 = fscanf(fid,'%f',1);
    disp(['Input Beta2:',num2str(Beta2)])
    if Beta2 <= 0
        disp('Error: Beta2 should > 0')
        return
    end
    if Beta1 >= Beta2
        disp('Error: Beta1 should < Beta2')
        return
    end
    Budget = fscanf(fid,'%f',1);
    disp(['Input Budget:',num2str(Budget)])
    if Budget <= 0
        disp('Error: Budget should > 0')
        return
    end
%     up = fscanf(fid,'%f',3);
    Epsilon = str2num(fscanf(fid,'%s',1));
    disp(['Input Epsilon:',num2str(Epsilon)])
    if Beta2 <= 0
        disp('Error: Epsilon should > 0')
        return
    end

    disp('Input Airport Security Check Conditions:')
    N = fscanf(fid,'%d',1);
    disp(['Input N:',num2str(N)])
    if N <= 0
        disp('Error: N should > 0')
        return
    end
    Lambda = str2num(fscanf(fid,'%s',1));
    disp(['Input Lambda:',num2str(Lambda)])
    if Lambda <= 0
        disp('Error: Lambda should > 0')
        return
    end
    
    mu1 = str2num(fscanf(fid,'%s',1));
    disp(['Input mu1:',num2str(mu1)])
    if mu1 <= 0
        disp('Error: mu1 should > 0')
        return
    end
    
    mu2 = str2num(fscanf(fid,'%s',1));
    disp(['Input mu2:',num2str(mu2)])
    if mu2 <= 0
        disp('Error: mu2 should > 0')
        return
    end
    if mu1 < mu2
        disp('Error: mu1 should > m2')
        return
    end
%     disp('Input theta:')
    
    theta = fscanf(fid,'%f',1);
    disp(['Input theta:',num2str(theta)])
    if theta <= 0
        disp('Error: theta should > 0')
        return
    end
	
% 	disp('Input d1')
    d1 = fscanf(fid,'%f',1);
    disp(['Input d1:',num2str(d1)])
    if d1 <= 0
        disp('Error: d1 should > 0')
        return
    end
	
% 	disp('Input d2:')
    d2 = fscanf(fid,'%f',1);
    disp(['Input d2:',num2str(d2)])
    if d2 <= 0
        disp('Error: d2 should > 0')
        return
    end
	
    disp('Select Method for Mean Waiting Time Evaluation')
    disp('1 ==> MatrixInverse, 2==> Chain')
    MeanWait_Flag = fscanf(fid,'%d',1);
    disp(['Input MeanWait_Flag:',num2str(MeanWait_Flag)])
    if MeanWait_Flag ~= 1 && MeanWait_Flag ~= 2
        disp('Selection should be 1 or 2')
        return
    end


    disp('Input SPLINE Parameters')
%     disp('Input Number of Initial Point')
    NInitial = fscanf(fid,'%d',1);
    disp(['Input Number of Initial Point:',num2str(NInitial)])
    if NInitial <= 0
        disp('Error: NInitial should > 0')
        return
    end

%     disp('Input the limit of traversed points')
    b = fscanf(fid,'%d',1);
    disp(['Input the limit of traversed points:',num2str(b)])
    if NInitial <= 0
        disp('Error: The limit of traversed points should > 0')
        return
    end
    
%     disp('Input the initial seed of initial points')
    seed = fscanf(fid,'%d',1);
    disp(['Input the initial seed of initial points:',num2str(seed)])
    if NInitial <= 0
        disp('Error: The seed should > 0')
        return
    end

    fclose(fid);
end