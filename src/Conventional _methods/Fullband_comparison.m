%% FULLBAND LINEAR COMPARISON METHODS

addpath '..\Common';                % Functions in Common folder
clear all;
close all;

%% Generate Response
M = 256;
mu = 0.001;
AR = 4;
SNR = 40;
iter = 1.0*80000;
b = load('h1.dat');              % Unknown system (select h1 or h2)
b = b(1:M);                      % Truncate to length M
% b = rand(1,M) - 0.5;

disp('Creating desired and input signals. . .');
fprintf('Input signal: noise colored wtih AR(%d) \n', AR);
fprintf('Linear system memory length: %d \n', M);

[un,dn,vn] = GenerateResponses(iter,b,sum(100*clock),AR,SNR);

% figures handlers
MSEfig = figure('Name', 'MSE');
Misfig = figure('Name', 'Misalingment');

fprintf('\n');

%% Wavelet Adaptive Filtering
fprintf('WAF \n');
fprintf('--------------------------------------------------------------------\n');
level = 3;                       % Levels of Wavelet decomposition
wtype = 'db4';                   % Mother Wavelet type
fprintf('Wavelet type: %s, levels: %d, step size = %.3f \n', wtype, level, mu);

S = WAFinit(zeros(M,1), mu, level, wtype);     % Initialization
S.unknownsys = b; 
tic;
[en, S] = WAFadapt(un, dn, S);                 % Perform WSAF Algorithm

fprintf('Total time = %.2f s \n',toc);

EML = S.eml.^2;                  % System error norm (normalized)
err_sqr = en.^2;
    

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

figure(MSEfig);                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE), 'DisplayName', 'WAF');

figure(Misfig);                          % Plot misalignment
hold on; plot((0:length(EML)-1)/1024,10*log10(EML),  'DisplayName', 'WAF');

fprintf('\n');

%% SOAF-DCT
fprintf('SOAF-DCT \n');
fprintf('--------------------------------------------------------------------\n');
fprintf('Step size: %.3f \n', mu);

S = SOAFinit(zeros(M,1),mu,iter);   % Initialization
S.unknownsys = b; 
tic;
[yn,en,S] = SOAFadapt(un,dn,S);     % Perform algorithm

fprintf('Total time = %.2f s \n',toc);

EML = S.eml.^2;                  % System error norm (normalized)
err_sqr = en.^2;
    

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

figure(MSEfig);                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE), 'DisplayName', 'SOAF-DCT');

figure(Misfig);                          % Plot misalignment
hold on; plot((0:length(EML)-1)/1024,10*log10(EML), 'DisplayName', 'SOAF-DCT');

fprintf('\n');

%% SOAF-DFT
% fprintf('SOAF-DFT \n');
% fprintf('--------------------------------------------------------------------\n');
% mu = 0.001;
% fprintf('Step size: %.4f \n', mu);
% 
% S = SOAFinit(zeros(M,1),mu,iter);   % Initialization
% S.unknownsys = b; 
% tic;
% [yn,en,S] = SOAFadapt_DFT(un,dn,S);     % Perform algorithm
% 
% fprintf('Total time = %.2f s \n',toc);
% 
% EML = abs(S.eml.^2);                  % System error norm (normalized)
% err_sqr = abs(en.^2);
%     
% 
% NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
% fprintf('NMSE = %.2f dB\n', NMSE);
% 
% figure(MSEfig);                          % Plot MSE
% q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
% hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE), 'DisplayName', 'SOAF-DFT');
% 
% figure(Misfig);                          % Plot misalignment
% hold on; plot((0:length(EML)-1)/1024,10*log10(EML), 'DisplayName', 'SOAF-DFT');
% 
% fprintf('\n');

%% Constrained FDAF
fprintf('Constrained FDAF \n');
fprintf('--------------------------------------------------------------------\n');
select = 1;
mu = 0.001;
mu_unconst = 0.001;
fprintf('Step size: %.3f \n', mu);

S = FDAFinit(zeros(M,1),mu,mu_unconst,iter);
S.select = select;
S.unknownsys = b;
tic;
[en,S] = FDAFadapt(un,dn,S);        % Perform FDAF algorithm

fprintf('Total time = %.2f s \n',toc);

EML = S.eml.^2;                  % System error norm (normalized)
err_sqr = en.^2;
   

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

figure(MSEfig);                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE), 'DisplayName', 'Constr. FDAF');

figure(Misfig);                          % Plot misalignment
hold on; plot((0:length(EML)-1)/1024,10*log10(EML), 'DisplayName', 'Constr. FDAF');

fprintf('\n');

%% Unconstrained FDAF
fprintf('Constrained FDAF \n');
fprintf('--------------------------------------------------------------------\n');
select = 2;
fprintf('Step size: %.3f \n', mu);

S = FDAFinit(zeros(M,1),mu,mu_unconst,iter);
S.select = select;
S.unknownsys = b;
tic;
[en,S] = FDAFadapt(un,dn,S);        % Perform FDAF algorithm

fprintf('Total time = %.2f s \n',toc);

EML = S.eml.^2;                  % System error norm (normalized)
err_sqr = en.^2;
   

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

figure(MSEfig);                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE), 'DisplayName', 'Unconstr. FDAF');

figure(Misfig);                          % Plot misalignment
hold on; plot((0:length(EML)-1)/1024,10*log10(EML), 'DisplayName', 'Unconstr. FDAF');

fprintf('\n');

%% NLMS
fprintf('NLMS \n');
fprintf('--------------------------------------------------------------------\n');
mu = 0.001;
fprintf('Step size: %.3f \n', mu);

S = NLMSinit(zeros(M,1),mu);        % Initialization
S.unknownsys = b;
tic;
[yn,en,S] = NLMSadapt(un,dn,S);     % Perform NLMS algorithm

fprintf('Total time = %.2f s \n',toc);

% EML = S.eml.^2;                  % System error norm (normalized)
err_sqr = en.^2;
   

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

figure(MSEfig);                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE), 'DisplayName', 'NLMS');

figure(Misfig);                          % Plot misalignment
hold on; plot((0:length(EML)-1)/1024,10*log10(EML), 'DisplayName', 'NLMS');

fprintf('\n');

%% Fixing plots...
figure(MSEfig);
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
title(sprintf('Input signal, colored with AR(%d)', AR));
legend('show');

figure(Misfig);  
axis([0 iter/1024 -50 5]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Misalignment (dB)'); title(sprintf('Input signal, colored with AR(%d)', AR));
grid on;
legend('show');
