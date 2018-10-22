% Testing Volterra filter
% clear all;
close all;

addpath "Fast Volterra Filtering";

fs = 1024;
N = 2*fs;
delta_w = 2/N;
w_freq = (0:N-1)/2; 
w_norm = (0:N-1)*delta_w;
t = 0:1/fs:2-1/fs;

freq = 10;
u = 1*sin(2*pi*freq*t);
% omega = pi/4*fs;
% u = 1*sin(omega*t);
U = fft(u, N);

M1 = 128; % length of first order volterra kernel
M2 = 64; % length of second order volterra kernel
gains = [1, 1];
ker1 = zeros(M1,1); 
ker2 = zeros(M2,M2); 


%Ker2 must be in triu form, if is tril only diagonal values are computed
ker1(1) = 1;
% ker1 = gains(2)*(rand(M1,1)-0.5);
ker2(1,1)= 1;
% ker2 = gains(2)*(rand(M2,M2)-0.5);
ker2 = eye(M2);
% ker2 = tril(ker2);
ker = {ker1, ker2};

y = fastVMcell(u, ker, [M1 M2]);
yr = sum(y,1);
yr = yr(1:length(u));
yr = yr - mean(yr);
Y = fft(yr, N);

%Plotting Stuffs
figure; 
subplot(3, 1, 1);
ylabel('Time (s)');
plot(t,u);

subplot(3, 1, 2);
plot(w_norm(1:N/2), abs(U(1:N/2)));
ylabel('Amplitude');

subplot(3, 1, 3);
plot(w_norm(1:N/2), abs(Y(1:N/2)));
xlabel('Normalized frequency \omega (\pi)'); ylabel('Amplitude');

figure; plot(yr);

