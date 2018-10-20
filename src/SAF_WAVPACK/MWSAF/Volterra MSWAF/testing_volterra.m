% Testing Volterra filter
clear all;
close all;

addpath "Fast Volterra Filtering";

fs = 1024;
N = 2*fs;
delta_w = 1/N;
w = (0:N-1);
freq = 100;
t = 0:1/fs:1-1/fs;
u = 1*sin(2*pi*freq*t);

figure; 
subplot(2, 2, [1 2]);
plot(t,u);

U = fft(u, N);

subplot(2, 2, [3 4]);
plot(w, abs(U));