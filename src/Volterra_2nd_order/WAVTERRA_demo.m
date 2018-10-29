% Volterra MWSAF       Multi Structured Wavelet-domain Adaptive Filter Demo
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]
diary log_WAVTERRA_TB.txt
fprintf('%s \n', datestr(datetime('now')));
addpath(genpath('../Common'));             % Functions in Common folder
addpath('DFT_bank Volterra'); 
addpath('../MWSAF'); 
clear all;  
close all;

%% Unidentified System parameters
order = 2; 
M1 = 256; % length of first order volterra kernel
M2 = 32; % length of second order volterra kernel

NL_system.M = [M1, M2];
gains = [1 0.1];

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
ker1 = rand(M1,1)-rand(1);
ker2 = second_order_kernel(M2);

%% Lowpass kernel
% norm_freq = 0.2;
% samples = [M1 M2]/2-1;
% b1 = norm_freq*sinc(norm_freq*(-samples(1)-1:samples(1)));
% b2 = norm_freq*sinc(norm_freq*(-samples(2)-1:samples(2)));
% ker1 = b1;
% ker2 = second_order_kernel(b2);

%% Simulated kernel - from h1 h2
% b1 = load('h1.dat');
% b1 = b1(1:M1);
% ker1 = b1;
% 
% b2 = load('h2.dat');
% b2 = b2(1:M2);
% ker2 = second_order_kernel(b2);


NL_system.Responses = {gains(1).*ker1, gains(2).*ker2};


% NL_system = create_volterra_sys(order, M, gains, 'nlsys1'); 

%% Plot 2-D kernel
kernel_plot(NL_system.Responses);

% Run parameters
iter = 1*80000;            % Number of iterations
disp('Creating desired and input signals. . .');
fprintf('Kernel Length: [%d, %d], iter= %d\n', M1, M2, iter);
[un,dn,vn] = GenerateResponses_Volterra(iter, NL_system ,sum(100*clock),1,40); %iter, b, seed, ARtype, SNR
% [un,dn,vn] = GenerateResponses_speech_Volterra(NL_system,'speech.mat');

figure;

for i = 1:2:10
%% Adaptive filter parameters
mu = [0.1, 0.1];            %Step sizes for different kernels 
C = M2;                     % Channels (kernel diagonals)

level = i;                  % Levels of Wavelet decomposition for different kernels
filters = 'db2';            % Set wavelet type for different kernels

SB = 1:2^level;                 % Nonlinear subband,

fprintf('Running iter %d of %d, level = %d , wtype= %s\n', i, 5, level, filters);

%% WAVTERRA
fprintf('--------------------------------------------------------------------\n');
fprintf('WAVTERRA\n');
fprintf('Wavelet type: %s, levels: %d, step size = %s \n', filters, level, sprintf('%s ', mu));

tic;
S = Volterra_Init(NL_system.M, mu, level, filters); 

% [en, S] = Volterra_2ord_adapt(un, dn, S);     
% [en, S] = Volterra_2ord_adapt_shift(un, dn, S, shift);   

S.true = NL_system.Responses; 
[en, S] = Volterra_2ord_adapt_v3(un, dn, S, C, SB);
% [en, S] = Volterra_2ord_adapt_opt(un, dn, S);

% [en, S] = Volterra_2ord_adapt_oldadapt(un, dn, S,10);

err_sqr = en.^2;
    
fprintf('Total time = %.2f s \n',toc);

% figure;         % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE), 'DisplayName', ['WAVTERRA - Level:', num2str(i), 'db4']);
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error'); grid on;
fprintf('NMSE = %.2f dB\n', 10*log10(sum(err_sqr)/sum(dn.^2)));
end

%%
fprintf('-------------------------------------------------------------\n');
fprintf('FULLBAND VOLTERRA NLMS\n');

% Run parameters
mu = [0.1, 0.1];

Sfull = Volterra_NLMS_init(NL_system.M, mu); 

tic;

[en, Sfull] = Volterra_NLMS_adapt(un, dn, Sfull);     

err_sqr_full = en.^2;
    
fprintf('Total time = %.2f s \n',toc);

% Plot MSE
q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr_full);
plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', 'FB NLMS');
axis([0 length(MSE_full)/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error'); grid on;
legend('show');

fprintf('NMSE = %.2f dB\n', 10*log10(sum(err_sqr_full)/sum(dn.^2)));

fprintf('\n');  % Empty line in logfile
diary off
