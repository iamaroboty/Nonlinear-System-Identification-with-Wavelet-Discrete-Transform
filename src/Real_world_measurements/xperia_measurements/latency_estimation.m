% latency estimation for xperia smartphone 
close all; 
clear all; 

imp = load('impulse.mat'); 
resp = load('resp.mat'); 

delta = imp.impulse; 
ir_res = resp.resp; 

plot(delta); 
hold on; 
plot(ir_res); 

fs = 44100 ; 


figure; 

N = fs; % number of FFT points

transform = fft(ir_res,N);
magTransform = 20*log10(abs(transform));

faxis = linspace(-fs/2,fs/2,N);
plot(faxis,fftshift(magTransform));
xlabel('Frequency (Hz)')

% view frequency content up to half the sampling rate:
axis([0 fs/2, -200 0]) 


%% latency 

[acorr, lag]= xcorr(delta, ir_res); 
[nn, i]= max(abs(acorr)); 

figure;
subplot(3,1,1); 
plot(delta);
subplot(3,1,2); 
plot(ir_res(1:size(delta,1)));
subplot(3,1,3); 
plot(lag, acorr); 



fprintf('latency is: %d \n', abs(lag(i)) ); 




