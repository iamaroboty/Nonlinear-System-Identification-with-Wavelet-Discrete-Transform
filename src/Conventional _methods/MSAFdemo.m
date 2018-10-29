% MSAFdemo         Multiband-Structured SAF (MSAF) Demo
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd

addpath '..\Common';             % Functions in Common folder
clear all;

% Adaptive filter parameters

mu = 0.1;                        % Step size (0<mu<2)
M = 256;                         % Length of adaptive weight vector
N = 4;                           % Number of subbands, 4
L = 8*N;                         % Length of analysis filters, M=2KN, 
                                 %   overlapping factor K=4

% Run parameters

iter = 1.0*80000;                % Number of iterations
b = load('h1.dat');              % Unknown system (select h1 or h2)
b = b(1:M);   
% b=zeros(M,1);
% b(1)=1;% Truncate to length M
% 
% % impulse response
% [y,Fs] = audioread('reverb_shimmer.wav');
% y = resample(y, 1,100);
% b = y(500:M+500-1)';

tic;

% Adaptation process

disp(sprintf('Number of subbands, N = %d, step size = %.2f',N,mu));
[un,dn,vn] = GenerateResponses(iter,b);
% [un,dn,vn] = GenerateResponses_speech(b,'SpeechSample.mat');
S = MSAFinit(zeros(M,1),mu,N,L); % Initialization
S.unknownsys = b;
[en,S] = MSAFadapt(un,dn,S);     % Perform MSAF algorithm

EML = S.eml.^2;                  % System error norm (normalized)
err_sqr = en.^2;
    
disp(sprintf('Total time = %.3f mins',toc/60));

figure;                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;

figure;                          % Plot misalignment
hold on; plot((0:length(EML)-1)/1024,10*log10(EML));
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Misalignment (dB)');
grid on;

