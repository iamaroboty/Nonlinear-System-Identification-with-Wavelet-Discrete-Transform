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
gains = [1 1];

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
kernel_plot(NL_system.Responses);


%% Adaptive filter parameters
mu = [0.1, 0.1];                 %Step sizes for different kernels 


level = [1];                  % Levels of Wavelet decomposition for different kernels
filters = 'db2';              % Set wavelet type for different kernels

level = 1;                  % Levels of Wavelet decomposition for different kernels
filters = 'db1';            % Set wavelet type for different kernels


% Run parameters
iter = 1.0*80000;                % Number of iterations

%%
tic;
% Adaptation process
fprintf('Wavelet type: %s, levels: %d, step size = %f \n', filters, level, sprintf('%f ', mu));
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
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE), 'DisplayName', 'Wavleterra');
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
fprintf('MSE = %.2f dB\n', mean(10*log10(MSE(end-2048:end))))


%% Fullband Volterra NLMS
fprintf('--------------------------------------------------------------------\n');
fprintf('FULLBAND VOLTERRA NLMS\n');
mu = [0.1, 0.1];
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


%% linear model
fprintf('--------------------------------------------------------------------\n');
fprintf('LINEAR MODEL\n');
mu = 0.05;
level = 2;
filters = 'db4';

Slin = SWAFinit(M1, mu, level, filters); 
[en, Slin] = MWSAFadapt(un, dn, Slin); 

err_sqr_lin = en.^2;
    
fprintf('Total time = %.3f mins \n',toc/60);

% Plot MSE
q = 0.99; MSE_lin = filter((1-q),[1 -q],err_sqr_lin);
hold on; plot((0:length(MSE_lin)-1)/1024,10*log10(MSE_lin), 'DisplayName', 'LWSAF');
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
fprintf('MSE_lin = %.2f dB\n', mean(10*log10(MSE_lin(end-2048:end))))
