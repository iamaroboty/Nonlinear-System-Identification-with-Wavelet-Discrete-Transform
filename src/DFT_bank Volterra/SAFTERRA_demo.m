% SAFdemo        Subband Adaptive Filter Demo
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd

addpath '..\Common';        % Functions in Common folder
clear all;


% Unidentified System parameters
order = 2; 
M1 = 256; % length of first order volterra kernel
M2 = 32; % length of second order volterra kernel

NL_system.M = [M1, M2];
gains = [1 1];

% Random Vector 
rng('default'); % For reproducibility

%% Simulated Kernel - random
ker1 = rand(M1,1) - rand(1);
ker2 = second_order_kernel(M2);

NL_system.Responses = {gains(1).*ker1, gains(2).*ker2};


% Adaptive filter parameters

mu = [0.1, 0.1];                   % Step size (0<mu<2)
M = [M1, M2];                    % Length of adaptive weight vector
N = 4;                      % Number of subbands, 4
D = N/2;                    % Decimation factor for 2x oversampling
L = 8*N;                    % Length of analysis filters, M=2KN, 
                            %   overlapping factor K=4

% Run parameters
iter = 1.0*80000;           % Number of iterations

%desired response gen
disp('Creating desired and input signals. . .');
[un,dn,vn] = GenerateResponses_Volterra(iter, NL_system ,sum(100*clock),1,40); %iter, b, seed, ARtype, SNR
%[un,dn,vn] = GenerateResponses_speech_Volterra(NL_system,'SpeechSample.mat');

%%
% Adaptation process
% Nonlinear model 
fprintf('--------------------------------------------------------------------\n');
fprintf('SAFTERRA\n');
disp(sprintf('Number of subbands, N = %d, step size = %.2f',N,mu));

S = SAFTERRA_Init(M,mu,N,D,L);
tic;
[en,S] = SAFTERRA_adapt(un,dn,S);
err_sqr = en.^2;

disp(sprintf('Total time = %.3f mins',toc/60));

figure;
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -80 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)');
grid on;

%%
% Nonlinear model 
fprintf('--------------------------------------------------------------------\n');
fprintf('MSAFTERRA\n');
disp(sprintf('Number of subbands, N = %d, step size = %.2f',N,mu));

S = MSAFTERRA_Init(M,mu,N,L);
tic;
[en,S] = MSAFTERRA_adapt(un,dn,S);
err_sqr = en.^2;

disp(sprintf('Total time = %.3f mins',toc/60));

q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -80 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)');
grid on;




