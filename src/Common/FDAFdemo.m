% FDAFdemo        Frequency-domain Adaptive Filter (FDAF) Demo
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd


addpath '..\Common';                % Functions in Common folder
clear all; close all;

% Adaptive filter parameters

M = 256;
mu = 0.1;                           % Step size for contrained FDAF 
mu_unconst = 0.001;                 % Step size for unconstrained FDAF 

% Run parameters

iter = 8.0*80000;                   % Number of iterations
un = randn(1,iter);                 % White noise

% Coloured noise

% un = randn(1,iter);               % Input signal un
% a = [1; -0.975;  0.95];           % AR(2) complex pole
% un = filter(1,a,un);
b = load('h1.dat');                 % Unknown system (select h1 or h2)
b = b(1:M);                         % Truncate to length M

select = input('Select (1) Constrained FDAF, (2) Unconstrained FDAF: ','s');
                                    % Select case (1) or (2) in run time

if str2num(select) == 1
    disp(sprintf('Constrained FDAF, step size = %.5f',mu));
else
    disp(sprintf('Unconstrained FDAF, step size = %.5f',mu_unconst));
end

% Adaptation process

[un,dn,vn] = GenerateResponses(iter,b);
tic;
S = FDAFinit(zeros(M,1),mu,mu_unconst,iter);
                                    % Initialization
S.unknownsys = b;
S.select = str2num(select);
[en,S] = FDAFadapt(un,dn,S);        % Perform FDAF algorithm

EML = S.eml.^2;                     % System error norm (normalized)
err_sqr = en.^2;
    
disp(sprintf('Total time = %.3f mins',toc/60));

figure;
temp = real(ifft(S.weight));
w = temp(1:S.length);
c = stem([b,w]);
legend('Actual','Estimated');
title('System identification of an FIR filter'); grid on;

figure;
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)');
title('FDAFdemo');
grid on;

figure;
hold on; plot((0:length(EML)-1)/1024,10*log10(EML));
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Misalignment (dB)');
title('FDAFdemo');
grid on;


