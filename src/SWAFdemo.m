% SWAFdemo          Subband Wavelet-domain Adaptive Filter Demo
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

addpath 'Common';             % Functions in Common folder
clear all;  close all;

% Adaptive filter parameters
mu = 0.4;                      % Step size
M = 256;                         % Length of unknown system response


level = 1;                       % Levels of Wavelet decomposition
wtype = 'db16';                   % Wavelet family
Ovr = 4;



% Run parameters
iter = 1.0*80000;                % Number of iterations
b = load('h1.dat');              % Unknown system (select h1 or h2)
b = b(1:M);                      % Truncate to length M


%b = sign(b);
% b = sign(b);



% TESTING, a = delay.
 %a = 2;
 %b = zeros(M,1);
 %b(a+1) = 1;

%% low pass filter system 
% norm_freq = 0.39;
% samples = M/2-1;
% 
% b = norm_freq*sinc(norm_freq*(-samples-1:samples));
%b = b + upsample(b(1:2:M).^4,2) + upsample(downsample(b,2).^6,2);
%b = horzcat(b, zeros(M-length(b)-1,1)');
% 
% %distort the low pass simple
% a = 0.5;
% k = 2*a/(1-a);
% b = (1+k)*(b)./(1+k*abs(b));

%  %% load reverb 
%  [y,Fs] = audioread('reverb_shimmer.wav');
%  M = length(y);
%  b = y(1:M);

%%
tic;
% Adaptation process
fprintf('Wavelet type: %s, levels: %d, step size = %f \n', wtype, level, mu);

[un,dn,vn] = GenerateResponses(iter,b,sum(100*clock),2,40); %iter, b, seed, ARtype, SNR
%S = SWAFinit(M, mu, level, wtype);   % Initialization
S = QMFInit(M, mu, level, wtype); 
S.unknownsys = b; 
[en, S] = SWAFadapt_DDDWT(un, dn, S);                 % Perform WSAF Algorithm 

err_sqr = en.^2;
    
fprintf('Total time = %.3f mins \n',toc/60);

figure;         % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
fprintf('MSE = %.2f dB\n', mean(10*log10(MSE(end-2048:end))))

%% time domain parameters
fs = 512; % samples per sec
freq = 100; % frequency
dt = 1/fs; 

%% impulse response
delta = [1; zeros(fs-2-1,1)];
figure; 
subplot(2,1,1)
stem(delta);
title('Input Signal'); 
% axis([0 10 -1.5 1.5])
out_resp = SWAFtest_DDDWT(delta, S); 
subplot(2,1,2)
stem(out_resp);
title('Output Signal-Estimated System vs True');
hold on; 
real_resp = filter(b, 1, delta);
stem(real_resp); 
% axis([0 2*M -1.5 1.5])

%% sine test signal 
amplitude = 1; 
leng = 1;
input_sine = amplitude*sin(2*pi*freq*(0:dt:leng-dt));

figure; 
subplot(2,2,1)
plot(input_sine);
title('Input Signal'); 
out_sine = SWAFtest_DDDWT(input_sine, S); 
subplot(2,2,2)
plot(out_sine);
title('Output Signal - Estimated System vs True');
hold on; 
real_sys = filter(b,1,input_sine);
plot(real_sys); legend('Estim', 'True');

%% FFT 
N = 2*fs;
faxis = linspace(-fs/2,fs/2,N);

subplot(2, 2, 3);
fft_true = abs(fft(input_sine, N)/N);
plot(faxis, fftshift(fft_true)); 
xlabel('Frequency');

subplot(2, 2, 4);
fft_out_est = abs(fft(out_sine, N)/N);
plot(faxis, fftshift(fft_out_est)); 
xlabel('Frequency');
hold on; 
fft_out_true = abs(fft(real_sys,N)/N);
plot(faxis, fftshift(fft_out_true));


% %%
% [C, L] = wavedec(b,level,wtype);
% cD = detcoef(C,L,level,'cell');
% cA = appcoef(C,L,wtype);
% 
% figure;
% % % Plot system ID differencies
% if level == 2
%     subplot(level +1, 1, 1);
%     stem([cD{1}, S.coeffs{1}]); legend('Actual','Estimated');
%     title('cD1'); grid on;
%     subplot(level +1, 1, 2);
%     stem([cD{2}, S.coeffs{2}(:,2)]); legend('Actual','Estimated');
%     title('cD2'); grid on;
%     subplot(level +1, 1, 3);
%     stem([cA, S.coeffs{2}(:,1)]); legend('Actual','Estimated');
%     title('cA'); grid on;
%     ax=axes('Units','Normal','Position',[.1 .1 .85 .85],'Visible','off');
%     set(get(ax,'Title'),'Visible','on')
%     title('System identification. Trasformed domain');
%     h=get(ax,'Title');
% 
%     b_est = waverec([S.coeffs{2}(:,1); S.coeffs{2}(:,2); S.coeffs{1}], S.L, wtype);
%     figure;
%     stem([b, b_est]);
%     legend('Actual','Estimated');
%     title('System identification. Time domain');grid on;
% 
% elseif level == 1
%     subplot(2, 1, 1);
%     stem([cD{1}, S.coeffs{1}(:,2)]); legend('Actual','Estimated');
%     title('cD'); grid on;
%     subplot(2, 1, 2);
%     stem([cA, S.coeffs{1}(:,1)]); legend('Actual','Estimated');
%     title('cA'); grid on;
%     ax=axes('Units','Normal','Position',[.1 .1 .85 .85],'Visible','off');
%     set(get(ax,'Title'),'Visible','on')
%     title('System identification. Trasformed domain');
%     h=get(ax,'Title');
% 
%     b_est = waverec([S.coeffs{1}(:,1); S.coeffs{1}(:,2)], S.L, wtype);
%     figure;
%     stem([b, b_est]);
%     legend('Actual','Estimated');
%     title('System identification. Time domain');grid on;
% else
%     fprintf('Set Level either 1 or 2\n');
% end
% 
% 
% 
