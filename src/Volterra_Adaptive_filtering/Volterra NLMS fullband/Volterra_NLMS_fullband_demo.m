% Volterra NLMS fullband demo       Fullband Volterra second order adaptive
%                                   filtering
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

diary log_VNLMS.txt
fprintf('%s \n', datestr(datetime('now')));
addpath(genpath('../../Common'));             % Functions in Common folder
clear all;  
close all;

%% Unidentified System parameters
order = 2; 
M1 = 256; % length of first order volterra kernel
M2 = 32; % length of second order volterra kernel

NL_system.M = [M1, M2];
gains = [1 2];

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
% N = 0; %diagonals number, beyond the main one
% for i = 1:N
%     d = diag(ones(M2-i,1),i);
%     ker2(d(:,:)==1) = rand(M2-i, 1) - rand(1);
% end    
% 
% d = eye(M2); ker2(d(:,:)==1) = rand(M2,1)- rand(1) ;     % instert principal diagonal

%% Simulated Kernel - random
% ker1 = rand(M1,1) - rand(1);
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

%% Fullband Volterra NLMS
fprintf('-------------------------------------------------------------\n');
fprintf('FULLBAND VOLTERRA NLMS\n');

% Run parameters
iter = 1.0*80000;                % Number of iterations
mu = [0.1, 0.1];

Sfull = Volterra_NLMS_init(NL_system.M, mu); 
[un,dn,vn] = GenerateResponses_Volterra(iter, NL_system ,sum(100*clock),4,40);
% [un,dn,vn] = GenerateResponses_speech_Volterra(NL_system,'speech.mat');

tic;

[en, Sfull] = Volterra_NLMS_adapt(un, dn, Sfull);     

err_sqr_full = en.^2;
    
fprintf('Total time = %.3f mins \n',toc/60);

% Plot MSE
figure;
q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr_full);
plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', 'FB NLMS');
axis([0 length(MSE_full)/1024 -inf 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
legend('show');

fprintf('NMSE = %.2f dB\n', 10*log10(sum(err_sqr_full)/sum(dn.^2)));

fprintf('\n');  % Empty line in logfile
diary off


