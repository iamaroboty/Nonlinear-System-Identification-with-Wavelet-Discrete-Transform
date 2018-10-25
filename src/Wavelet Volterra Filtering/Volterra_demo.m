% Volterra MWSAF       Multi Structured Wavelet-domain Adaptive Filter Demo
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

addpath '../../../../Common';             % Functions in Common folder
clear all;  
close all;

%% Unidentified System parameters
order = 2; 
M1 = 256; % length of first order volterra kernel
M2 = 32; % length of second order volterra kernel

NL_system.M = [M1, M2];
gains = [1 0.5];

%NL_system = create_volterra_sys(order, M, gains, 'nlsys1'); 
%% Just a Delta
% ker1 = zeros(M1,1);
% ker1(1) = 1;
% % ker2 = diag(ones(M2,1));
% ker2 = zeros(M2,M2);
% ker2(1,1) = 2;

%% Random Vector 
rng('default'); % For reproducibility
% ker1 = rand(M1,1) - rand(1);
% shift = 0;      
% ker2 = diag(rand(M2-shift,1)- rand(1) , shift); 
% 
% N = 5; %diagonals number, beyond the main one
% for i = 1:N
%     d = diag(ones(M2-i,1),i);
%     ker2(d(:,:)==1) = rand(M2-i, 1) - rand(1);
% end    
% 
% d = eye(M2); ker2(d(:,:)==1) = rand(M2,1)- rand(1) ;     % instert principal diagonal

%% Simulated Kernel - random
% ker1 = rand(M1,1) - 0.5;
% ker2 = second_order_kernel(M2);

%% Simulated kernel - from h1 h2
b1 = load('h1.dat');
b1 = b1(1:M1);
ker1 = b1;

b2 = load('h2.dat');
b2 = b2(1:M2);
ker2 = second_order_kernel(b2);



NL_system.Responses = {gains(1).*ker1, gains(2).*ker2};

% NL_system = create_volterra_sys(order, M, gains, 'nlsys1'); 

%% Plot 2-D kernel
figure; 
subplot(221);   %linear time
stem(NL_system.Responses{1});
ylabel('Linear Kernel (Time)');
xlabel('Sample (N)');
grid on;

subplot(222);   % linear freq
n_points = 1024;
w = linspace(-1,1, n_points);
K1 = fft(NL_system.Responses{1}, n_points);
plot(w, 20*log10(fftshift(abs(K1))));
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Linear Kernel - Power (dB)');
grid on;

sub3 = subplot(223);   % quadratic time
surf(NL_system.Responses{2});
ylabel('Sample (N)');
xlabel('Sample (N)');
view([0 90]);
colormap(sub3, parula);
colorbar
grid on;

sub4 = subplot(224);   % quadratic freq
title('Quadratic Kernel (Frequency)');
K2 = fft2(NL_system.Responses{2}, n_points, n_points);
surf(w, w, 20*log10(fftshift((abs(K2)))), 'LineStyle', 'none');
zlabel('Power (dB)');
ylabel('Normalized Frequency (\times\pi rad/sample)');
xlabel('Normalized Frequency (\times\pi rad/sample)');
view([0 90]);
colormap(sub4, jet);
colorbar
grid on;


%% Adaptive filter parameters
mu = [0.1, 0.1];                 % Step sizes for different kernels 

level = [4];                  % Levels of Wavelet decomposition for different kernels
filters = 'db8';               % Set wavelet type for different kernels

% Run parameters
iter = 2.0*80000;                % Number of iterations

%%
tic;
% Adaptation process
fprintf('Wavelet type: %s, levels: %d, step size = %f \n', filters, level, mu);
[un,dn,vn] = GenerateResponses_Volterra(iter, NL_system ,sum(100*clock),1,40); %iter, b, seed, ARtype, SNR
% [un,dn,vn] = GenerateResponses_speech_Volterra(NL_system,'SpeechSample.mat');

%% Nonlinear model 
S = Volterra_Init(NL_system.M, mu, level, filters); 

% [en, S] = Volterra_2ord_adapt(un, dn, S);     
% [en, S] = Volterra_2ord_adapt_shift(un, dn, S, shift);   
[en, S] = Volterra_2ord_adapt_v2(un, dn, S);


err_sqr = en.^2;
    
fprintf('Total time = %.3f mins \n',toc/60);

figure;         % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE), 'DisplayName', 'WavTerra');
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
fprintf('MSE = %.2f dB\n', mean(10*log10(MSE(end-2048:end))))


%% Fullband Volterra NLMS
fprintf('-------------------------------------------------------------\n');
fprintf('FULLBAND VOLTERRA NLMS\n');
mu = [0.3, 0.1];
Sfull = Volterra_NLMS_init(NL_system.M, mu); 

[en, Sfull] = Volterra_NLMS_adapt(un, dn, Sfull);     

err_sqr_full = en.^2;
    
fprintf('Total time = %.3f mins \n',toc/60);

% Plot MSE
q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr_full);
hold on; plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', 'FB NLMS');
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
fprintf('MSE = %.2f dB\n', mean(10*log10(MSE_full(end-2048:end))))
legend('show');


% %% linear model
% 
% fprintf('-----------------------------------------------------------\n');
% fprintf('LINEAR MODEL\n');
% 
% 
% S = SWAFinit(M1, mu(1), level(1), filters); 
% [en, S] = MWSAFadapt(un, dn, S); 
% 
% err_sqr = en.^2;
%     
% fprintf('Total time = %.3f mins \n',toc/60);
% 
% figure;         % Plot MSE
% q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
% hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
% axis([0 iter/1024 -60 10]);
% xlabel('Number of iterations (\times 1024 input samples)'); 
% ylabel('Mean-square error (with delay)'); grid on;
% fprintf('MSE = %.2f dB\n', mean(10*log10(MSE(end-2048:end))))



% figure;                          % Plot misalignment
% hold on; plot((0:length(EML)-1)/1024,10*log10(EML));
% xlabel('Number of iterations (\times 1024 input samples)'); 
% ylabel('Misalignment (dB)');
% grid on;
