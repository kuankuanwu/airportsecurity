theta=0.35;

alpha_vector=0:0.005:1;
f_vector = zeros(length(alpha_vector),1); 

for i = 1 : length(alpha_vector)
    f_vector(i) = exp(-alpha_vector(i)/theta)/(theta*(1-exp(-1/theta))); 
end

plot(alpha_vector,f_vector);

tau_vector = 0:0.005:1; 
p_vector = zeros(length(tau_vector),1);

for i = 1 : length(tau_vector)
    p_vector(i) = (-exp(-tau_vector(i)/theta)+1)/(1-exp(-1/theta));
end

plot(tau_vector,p_vector);