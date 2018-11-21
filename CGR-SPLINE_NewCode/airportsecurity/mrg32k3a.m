function [rn, s22] = mrg32k3a(s22)

 norm = 2.328306549295728E-10;
 m1   = 4294967087;
  m2   = 4294944443;
  a12  = 1403580;
             a13n = 810728;
             a21  = 527612;
             a23n = 1370589;

 s10=1;
 s11=2;
 s12=3;
 s20=4;
 s21=5;
             
if s22 <=0
   
    s22 = -s22;
    s10 = s22+1;
    
    if s10>=m1 
       s10 = 1; 
    end    
    
    s11=s22+2;
    
    if s11>=m1
       s11 = 2;
        
    end    
    
    s12=s22+3;
    
    if s12>=m1
        s12=3;
    end
    
    s20=s22+4;
    
    if s20>=m2
        
        s20 = 4;
    end 
    
    s21=s22+5;
    if s21>=m2
        
        s21=5;
    end
    
end

p1=a12*s11-a13n*s10;
k = fix(p1/m1);
p1=p1-k*m1;

if p1<0
    p1 = p1+m1;
    
end

s10=s11;
s11=s12;
s12=p1;

p2=a21*s22-a23n*s20;
k=fix(p2/m2);
p2=p2-k*m2;

if p2<=0
    
    p2=p2+m2;
end

s20=s21;
s21=s22;
s22=p2;

if p1<=p2
    
    rn = (p1-p2+m1)*norm;
    
else
    rn = (p1-p2)*norm;
end

end