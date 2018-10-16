% SAFdemo        Subband Adaptive Filter Demo
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd

addpath '..\Common';        % Functions in Common folder
clear all;

% Adaptive filter parameters

mu = 0.1;                   % Step size (0<mu<2)
M = 256;                    % Length of adaptive weight vector
N = 4;                      % Number of subbands, 4
D = N/2;                    % Decimation factor for 2x oversampling
L = 8*N;                    % Length of analysis filters, M=2KN, 
                            %   overlapping factor K=4

% Run parameters

iter = 1.0*80000;           % Number of iterations
b = load('h1.dat');         % Unknown system (select h1 or h2)
b = b(1:M);                 % Truncate to length M

%b = sign(b)
%b = tanh(b.*100);
%b=zeros(1,M);
%b(64)=1;


% Adaptation process

disp(sprintf('Number of subbands, N = %d, step size = %.2f',N,mu));
%[un,dn,vn] = GenerateResponses(iter,b);
[un,dn,vn] = GenerateResponses_speech(b, 'SpeechSample.mat');
S = SAFinit(M,mu,N,D,L);
tic;
[en,S] = SAFadapt(un,dn,S);
err_sqr = en.^2;

disp(sprintf('Total time = %.3f mins',toc/60));

figure;
q = 0.9999; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -80 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)');
grid on;


