clc

InputFileName = 'input.txt';
InitializeAndReadInput(InputFileName)

S1=2; 
S2=3;

tau_vector = 0:0.005:1;
meanwait_vector = zeros(length(tau_vector),1);

for i = 1 : length(tau_vector)
    meanwait_vector(i) = MeanWait(tau_vector(i),S1,S2);
end

plot(tau_vector,meanwait_vector); 