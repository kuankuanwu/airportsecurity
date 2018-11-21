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
    global b
    global d1
    global d2

    global fout
    fid = fopen(InputFileName,'r');
    fout = fopen('Output.txt','w');

    %(Not Critical Parameters for Finding Minimal tau)
    d1 = 0.7;
    d2 = 0.98;


    if fid < 0
        disp('Error: File Not Found!');
        return
    end
    disp('Input Constraint Conditions:')
    disp('Input c1:')
    c1 = fscanf(fid,'%f',1);
    if c1 <= 0
        disp('Error: c1 should > 0')
        return
    end
    disp('Input c2:')
    c2 = fscanf(fid,'%f',1);
    if c2 <= 0
        disp('Error: c2 should > 0')
        return
    end
    if c1 >= c2
        disp('Error: c1 should < c2')
        return
    end
    disp('Input Beta1:')
    Beta1 = fscanf(fid,'%f',1);
    if Beta1 <= 0
        disp('Error: Beta1 should > 0')
        return
    end
    disp('Input Beta1:')
    Beta2 = fscanf(fid,'%f',1);
    if Beta2 <= 0
        disp('Error: Beta2 should > 0')
        return
    end
    if Beta1 >= Beta2
        disp('Error: Beta1 should < Beta2')
        return
    end
    disp('Input Budget:')
    Budget = fscanf(fid,'%f',1);
    if Budget <= 0
        disp('Error: Budget should > 0')
        return
    end
    disp('Input Epsilon:')
    Epsilon = fscanf(fid,'%f',1);
    if Beta2 <= 0
        disp('Error: Epsilon should > 0')
        return
    end

    disp('Input Airport Security Check Conditions:')
    disp('Input N:')
    N = fscanf(fid,'%d',1);
    if N <= 0
        disp('Error: N should > 0')
        return
    end
    disp('Input Lambda:')
    Lambda = fscanf(fid,'%f',1);
    if Lambda <= 0
        disp('Error: Lambda should > 0')
        return
    end
    disp('Input mu1:')
    mu1 = fscanf(fid,'%f',1);
    if mu1 <= 0
        disp('Error: mu1 should > 0')
        return
    end
    disp('Input mu2:')
    mu2 = fscanf(fid,'%f',1);
    if mu2 <= 0
        disp('Error: mu2 should > 0')
        return
    end
    if mu1 < mu2
        disp('Error: mu1 should > m2')
        return
    end
    disp('Input theta:')
    theta = fscanf(fid,'%f',1);
    if theta <= 0
        disp('Error: theta should > 0')
        return
    end

    disp('Select Method for Mean Waiting Time Evaluation')
    disp('1 ==> MatrixInverse, 2==> Chain')
    MeanWait_Flag = fscanf(fid,'%d',1);
    if MeanWait_Flag ~= 1 && MeanWait_Flag ~= 2
        disp('Selection should be 1 or 2')
        return
    end


    disp('Input SPLINE Parameters')
    disp('Input Number of Initial Point')
    NInitial = fscanf(fid,'%d',1);
    if NInitial <= 0
        disp('Error: NInitial should > 0')
        return
    end

    disp('Input the limit of traversed points')
    b = fscanf(fid,'%d',1);
    if NInitial <= 0
        disp('Error: The limit of traversed points should > 0')
        return
    end

    fclose(fid);
end