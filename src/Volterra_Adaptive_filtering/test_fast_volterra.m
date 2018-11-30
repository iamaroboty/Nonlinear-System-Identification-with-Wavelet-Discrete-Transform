%% fast volterra filtering test
fs = 8000; 
freq = 20; 
dt = 1/fs; 

amplitude = 1; 
leng = 1; 

s = amplitude*sin(2*pi*freq*(0:dt:leng-dt));
s2 = s + s.^2;

plot(s2); 
hold on; 

ker1 = zeros(256,1); 
ker1(1) = 1;
ker2 = zeros(16,16); 
ker2(1,1) =1; 
kernels = {ker1,ker2};
y = sum(fastVMcell(s, kernels, [256,16] ),1); 

plot(y)

