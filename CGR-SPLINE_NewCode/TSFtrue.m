function TSFtrue = TSFtrue(problemparam, x)

% =========================================================================
% Calculates throughput *only* (analytically) for the the three-stage flowline 
% problem. (Reference ?)
%
% problemparam = 
%    param(1) = problem ID, not used
%    param(2) = problem dimension = 4
%    param(3) = nseeds = 3
%    param(4) = nSecMeas = 1 
%    param(5) = warm-up time
%    param(6) = simulation end time
%    param(7) = total service rate
%    param(8) = total buffer space available
% x =
%	 x(0) = id
%	 x(1) = mu1
%	 x(2) = mu2
%	 x(3) = mu3
%    x(4) = b2
%
% Example: 
% param = [123 4 3 1 50 1000 20 20]
% x=[6 7 8 9]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ratetotal = int32(problemparam(7));
buffertotal=int32(problemparam(8));
mu1 = x(1);
mu2 = x(2);
mu3 = x(3); %//rate - mu1 - mu2;

if sum(x<=0)>0 || sum(x(1:3)>ratetotal)>0 || x(4)>buffertotal-1 
    TSFtrue=-9999;
    return
end

b2 = x(4); %//x(3);
b3 = buffertotal - b2;
id2 = (b2+2)*(b3+2) - 1;
	
	%fprintf(1, 'b2 = %d, b3 = %d\n', b2, b3);
	%fprintf(1, 'mu1 = %d, mu2 = %d, mu3= %d\n', mu1, mu2, mu3);
	%fprintf(1, 'id2 = %d\n', id2);

a=zeros(id2, id2);
b=zeros(id2, 1);

%//for p(i,j), column = i*(b3+2) + j
%//row = 0, column = i*(b3+1) + j (for p(i,j))
%//mu1*p(0,0) - mu2*p(0,1) = 0
a(0+1, 0+1) = mu1;
a(0+1, 1+1) = -mu3;
b(0+1)      = 0;

%//row = 1:b3
%//(mu1 + mu3)*p(0,n3) - mu2*p(1, n3-1) - mu3*p(0,n3+1) = 0; 1 <= n3 <= b3 
for k=1:b3
	a(k+1, k+1)             = mu1 + mu3;
	a(k+1, b3 + k + 1 + 1)  = -mu2;
	a(k+1, k+1+1)           = -mu3;
	b(k+1)                  = 0;
end
	
%//row = b3+1
%//(mu1 + mu3)*p(0,b3+1) - mu2*p(1, b3) = 0 	
a(b3+1+1, b3+1+1)   = mu1 + mu3;
a(b3+1+1, 2*b3+2+1) = -mu2;
b(b3+1+1)           = 0;
	

%//row = (b3+1)+1:(b3+1)+b2 **
%//(mu1+mu2)*p(n2, 0) - mu1*p(n2-1,0) - mu3*p(n2,1); 1 <= n2 <= b2
for k=1:b2
	a(b3+1+k+1, k*(b3+2)+1)       = mu1 + mu2;
	a(b3+1+k+1, (k-1)*(b3+2)+1)   = -mu1;
	a(b3+1+k+1, k*(b3+2)+1+1)     = -mu3;
	b(b3+1+k+1)                   = 0;
end


%//row = (b2+b3+1)+1:(b2+b3+1)+b2*b3
%//(mu1 + mu2 + mu3)*p(n2, n3) - mu1*p(n2-1, n3) - mu2*p(n2+1, n3-1) - mu3*p(n2, n3+1) = 0; 1 <= n2 <= b2, 1 <= n3 <= b3
for i=1:b2
    for j=1:b3
		a(b2+b3+1+(i-1)*b3+j+1, i*(b3+2)+j+1)          = mu1 + mu2 + mu3;
		a(b2+b3+1+(i-1)*b3+j+1, (i-1)*(b3+2)+j+1)      = -mu1;
		a(b2+b3+1+(i-1)*b3+j+1, (i+1)*(b3+2)+(j-1)+1)  = -mu2;
		a(b2+b3+1+(i-1)*b3+j+1, i*(b3+2)+(j+1)+1)      = -mu3;
		b(b2+b3+1+(i-1)*b3+j+1)                        = 0;
    end
end
	
		
%//row = (b2*b3+b2+b3+1)+1:(b2*b3+b2+b3+1)+b2-1
%//(mu1+mu3)*p(n2,b3+1) - mu1*p(n2-1,b3+1) - mu2*p(n2+1,b3) = 0; 1 <= n2 <= b2-1
for k=1:b2-1
	a(b2*b3+b2+b3+1+k+1, k*(b3+2)+b3+1+1)       = mu1 + mu3;
	a(b2*b3+b2+b3+1+k+1, (k-1)*(b3+2)+b3+1+1)   = -mu1;
	a(b2*b3+b2+b3+1+k+1, (k+1)*(b3+2)+b3+1)     = -mu2;
	b(b2*b3+b2+b3+1+k+1)                        = 0;
end
 

%//row = (b2*b3+2*b2+b3)+1
%//mu3*p(b2,b3+1) - mu1*p(b2-1,b3+1) - mu2*p(b2+1,b3)
a((b2*b3+2*b2+b3)+1+1, b2*(b3+2)+(b3+1)+1)         = mu3;
a((b2*b3+2*b2+b3)+1+1, ((b2-1)*(b3+2))+(b3+1)+1)   = -mu1;
a((b2*b3+2*b2+b3)+1+1, (b2+1)*(b3+2)+b3+1)         = -mu2;
b((b2*b3+2*b2+b3)+1+1)                             = 0;
	

%//row = (b2*b3+2*b2+b3+1)+1
%//mu2*p(b2+1,0) - mu1*p(b2,0) - mu3*p(b2+1,1) = 0
a((b2*b3+2*b2+b3+1)+1+1, (b2+1)*(b3+2)+1)    = mu2;
a((b2*b3+2*b2+b3+1)+1+1, b2*(b3+2)+1)        = -mu1;
a((b2*b3+2*b2+b3+1)+1+1, (b2+1)*(b3+2)+1+1)  = -mu3;
b((b2*b3+2*b2+b3+1)+1+1)                     = 0;


%//row = (b2*b3+2*b2+b3+2)+1:(b2*b3+2*b2+b3+2)+b3-1
%//(mu2+mu3)*p(b2+1,n3) - mu1*p(b2,n3) - mu3*p(b2+1,n3+1) = 0; 1 <= n3 <= b3-1
for k=1:b3-1
	a(b2*b3+2*b2+b3+2+k+1, (b2+1)*(b3+2)+k+1)       = mu2 + mu3;
	a(b2*b3+2*b2+b3+2+k+1, b2*(b3+2)+k+1)           = -mu1;
	a(b2*b3+2*b2+b3+2+k+1, (b2+1)*(b3+2)+(k+1)+1)   = -mu3;
	b(b2*b3+2*b2+b3+2+k+1)                          = 0;
end


%//row = b2*b3+2*(b2+b3)+2
%//sum_{(n2,n3) in S} p(n2,n3) = 1; 
for k=0:(b3+2)*(b2+2)-2
    a(b2*b3+2*(b2+b3)+2+1, k+1) = 1;
end	
b(b2*b3+2*(b2+b3)+2+1) = 1;

x = linsolve(a, b);

    %fprintf(1, 't = ( ');
    %for i=0:id2-1
    %    fprintf(1, '%.6f ', x(i+1));
    %end
    %fprintf(1, ')\n');


TSFtrue=0;
for i=0:b2+1
	TSFtrue = TSFtrue + x(i*(b3+2)+1);
	%printf("%d %.20f\n", i*(b3+2), t.var(i*(b3+2)));
end
TSFtrue = mu3*(1-TSFtrue);
	%fprintf(1, 'Last constraint: %.20f\n', (mu2+mu3)*x((b2+1)*(b3+2)+b3+1)-mu1*x(b2*(b3+2)+b3+1));
    %fprintf(1, 'throughput=%.12f\n', thruput);