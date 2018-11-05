% SOAFdemo        Self-orthogonalizing Adaptive Filter (SOAF) Demo using DCT 
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd


addpath '..\Common';                % Functions in Common folder
clear all; close all;

% Adaptive filter parameters

M = 256;                            % Length of adaptive filter
mu = 0.001;                         % Default step size

% Run parameters

iter = 8.0*80000;                   % Number of iterations at each run
%un = randn(1,iter);                % White noise
%Coloured noise
un = randn(1,iter);                 % Generate random signal
a = [1; -0.975;  0.95];             % AR(2) model (complex poles)
un = filter(1,a,un);                % Generate AR signal
b = load('h1.dat');                 % Unknown system (select h1 or h2)
b = b(1:M);                         % Truncate to length M



% Adaptation process

disp(sprintf('SOAF, step size = %.5f',mu));
[un,dn,vn] = GenerateResponses(iter,b);
tic;
S = SOAFinit(zeros(M,1),mu,iter);   % Initialization
S.unknownsys = b;
[yn,en,S] = SOAFadapt(un,dn,S);     % Perform algorithm

EML = S.eml.^2;                     % System error norm (normalized)
err_sqr = en.^2;
    
disp(sprintf('Total time = %.3f mins',toc/60));

figure;
B = S.T*b;
c = stem([B,S.coeffs]);
legend('Actual','Estimated');
title('System identification of FIR filter');grid on;

figure;                             % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (dB)');
title('SOAFdemo');
grid on;

figure;                             % Plot misalignment
hold on; plot((0:length(EML)-1)/1024,10*log10(EML));
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Misalignment (dB)');
title('SOAFdemo');
grid on;


