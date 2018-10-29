% NLMSdemo         Normalized LMS (NLMS) Algorithm Demo
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd

addpath '..\Common';                % Functions in Common folder
clear all;

% Adaptive filter parameters

mu = 0.1;                           % Default step size (0<mu<2)
M = 256;                            % Length of adaptive weight vector

% Run parameters

iter = 1.0*80000;                   % Number of iterations
b = load('h1.dat');                 % Unknown system (select h1 or h2)
b = b(1:M);                         % Truncate to length M


% Adaptation process

disp(sprintf('NLMS, step size = %.2f',mu));

%[un,dn,vn] = GenerateResponses(iter,b);

[un,dn,vn] = GenerateResponses_speech(b, 'SpeechSample.mat');
tic;
S = NLMSinit(zeros(M,1),mu);        % Initialization
S.unknownsys = b;
[yn,en,S] = NLMSadapt(un,dn,S);     % Perform NLMS algorithm

EML = S.eml.^2;                     % System error norm (normalized)
err_sqr = en.^2;
    
disp(sprintf('Total time = %.3f mins',toc/60));

figure;
q = 0.9999; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -80 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (dB)');
title('NLMSdemo');
grid on;

figure;
hold on; plot((0:length(EML)-1)/1024,10*log10(EML));
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Misalignment (dB)');
title('NLMSdemo');
grid on;

