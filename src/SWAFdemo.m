% SWAFdemo          Subband Wavelet-domain Adaptive Filter Demo
% 
% by A. Castellani & S. Cornell [Università Politecnica delle Marche]

addpath 'Common';             % Functions in Common folder
clear all; close all;

% Adaptive filter parameters

mu = 0.1;                      % Step size
M = 256;                         % Length of unknown system response
level = 2;                       % Levels of Wavelet decomposition
wtype = 'db1';                   % Wavelet family

% Run parameters
iter = 1.0*80000;                % Number of iterations
b = load('h1.dat');              % Unknown system (select h1 or h2)
b = b(1:M);                      % Truncate to length M
% b =[1; zeros(M-2,1)];
% b =[zeros(M/2-1,1); 1; zeros(M/2,1)];
tic;

% Adaptation process

fprintf('Wavelet type: %s, levels: %d, step size = %f \n', wtype, level, mu);
[un,dn,vn] = GenerateResponses(iter,b,sum(100*clock),1,40); %iter, b, seed, ARtype, SNR
S = SWAFinit(M, mu, level, wtype);     % Initialization
S.unknownsys = b; 
[en, S] = SWAFadapt(un, dn, S);                 % Perform WSAF Algorithm

err_sqr = en.^2;
    
fprintf('Total time = %.3f mins \n',toc/60);

figure;                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
