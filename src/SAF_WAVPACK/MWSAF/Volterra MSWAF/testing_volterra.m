% Testing Volterra filter
clear all;
close all;

addpath "Fast Volterra Filtering";

fs = 1024;
N = 2*fs;
delta_w = 2/N;
w_freq = (0:N-1)/2; 
w_norm = (0:N-1)*delta_w;
freq = 10;
t = 0:1/fs:1-1/fs;

u = 1*sin(2*pi*freq*t);

figure; 
subplot(3, 1, 1);
plot(t,u);

U = fft(u, N);

subplot(3, 1, 2);
plot(w_freq(1:N/2), abs(U(1:N/2)));
xlabel('Normalized frequency \omega (\pi)'); ylabel('Amplitude');

M1 = 256; % length of first order volterra kernel
M2 = 32; % length of second order volterra kernel
gains = [1, 1];
ker1 = zeros(M1,1); 
ker1(1) = gains(2);
ker2 = zeros(M2,M2); 
ker2(1,1) = gains(2);
ker = {ker1, ker2};
y = fastVMcell(u, ker, [M1 M2]);
yr = sum(y,1);

subplot(3, 1, 3);
plot(abs(fft(yr, N)));

figure; plot(abs(fft(y(2,:))));

