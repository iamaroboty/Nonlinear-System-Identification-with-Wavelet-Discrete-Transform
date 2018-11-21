% Linear model comparison
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

addpath '../Common';             % Functions in Common folder

clear all;  
close all;
rng('default');

%% Structure Hyperparameters
% Filters length
M = 512;

% Universal stepsize
mu = 0.3;
mu_small = 0.003;

% For Wavelet Transform
level = 3;
filters = 'db8';

% Hammerstein Full Band
% p , w
order = 1;
mu_p = 0.3;
mu_w = 0.5;
alpha = 10.^0;

% For MSAF and SAF
N = 2^level;                % Number of subbands                     
D = N/2;                    % Decimation factor for 2x oversampling
L = 8*N;                    % Length of analysis filters, M=2KN


% Iter count
iter = 1.0*80000;

infile = 'color';
gain = '5p';

%% Input Data
% Anecoic measurment
load(['ref_', infile, '.mat'])
un = un(:)';    % Make it row vector
load(['univpm_', infile, '_', gain,'.mat'])
dn = dn(:)';    % Make it row vector

delay = 1024 + 100;
dn = dn(delay:end); 
un = un(1:end-delay+1); 

dn = dn(1:iter);
un = un(1:iter);

fs = 44100;

% % Resample the data to new fsn
% fs_new = 16000;
% [P, Q] = rat(fs_new/fs);
% dn = resample(dn, P, Q);
% un = resample(un, P, Q);

% Normalization of u(n) and d(n)
un = un/std(dn);
dn = dn/std(dn);

% % Normalization of speech signal
% a = abs(max(dn))/sqrt(2);
% un = un/a; dn = dn/a;

%%
% [un,dn,vn] = GenerateResponses(iter,b,sum(100*clock),4,40); %iter, b, seed, ARtype, SNR
% [un,dn,vn] = GenerateResponses_speech(b,'SpeechSample.mat');


%%
iter = length(un);
fprintf('Prior assumed filter Lengths: %d, iter= %d\n', M, iter);

% figures handlers
MSEfig = figure('Name', 'MSE');
NMSEfig = figure('Name', 'NMSE');

fprintf('\n');

%% NLMS
fprintf('NLMS \n');
fprintf('--------------------------------------------------------------------\n');
fprintf('Step size: %.3f \n', mu_small);

S = NLMSinit(zeros(M,1),mu_small);        % Initialization
tic;
[~,en,S] = NLMSadapt(un,dn,S);     % Perform NLMS algorithm

fprintf('Total time = %.2f s \n',toc);

err_sqr = en.^2;
   
NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

figure(MSEfig);                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:iter-1)/1024,10*log10(MSE), 'DisplayName', 'NLMS');

figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'NLMS' );


fprintf('\n');

%% WMSAF
fprintf('WMSAF \n');
fprintf('--------------------------------------------------------------------\n');
fprintf('Wavelet type: %s, levels: %d, step size = %f \n', filters, level, mu);

tic;
S = SWAFinit(M, mu, level, filters); 
[en, S] = MWSAFadapt(un, dn, S);    
%     [en_cllp, S_cllp] = MWSAFadapt_cllp(un, dn, S);  
%     [en_oplp, S_oplp] = MWSAFadapt_oplp(un, dn, S);  

fprintf('Total time = %.2f s \n',toc);

err_sqr = en.^2;
   
NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

figure(MSEfig);                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:iter-1)/1024,10*log10(MSE), 'DisplayName', 'WMSAF');

figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'WMSAF' );


fprintf('\n');

%% HAMMERESTEIN
fprintf('HAMMERSTEIN\n');
fprintf('--------------------------------------------------------------------\n');
fprintf('Order = %d, sys_len = %d, mu = %d \n', order, M, mu);        

tic;
SH = Hammerstein_NLMS_init(order, M, [mu_p mu_w], [0 0], alpha); 
[en, SH] = Hammerstein_NLMS_adapt(un, dn, SH);     

err_sqr = en.^2;

fprintf('Total time = %.2f s \n',toc);

figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:iter-1)/1024,10*log10(MSE), 'DisplayName', 'HAMM Or.1');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'HAMM Or.1' );

fprintf('\n');

%% WAMMERSTERIN
fprintf('WAMMERSTERIN\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Order = %d, sys_len = %d, mu = %d \n', order, M, mu);        

tic;
S = Wammerstein_init(M, [mu_p mu_w] ,level, filters, order, alpha); 

[en, S] = Wammerstein_adapt(un, dn, S);  

err_sqr = en.^2;

fprintf('Run time = %.2f s \n',toc);

figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:iter-1)/1024,10*log10(MSE), 'DisplayName', 'WAMM Or.1');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'WAMM Or.1' );

fprintf('\n');

%% SAF
fprintf('SAF\n');
fprintf('-------------------------------------------------------------\n');

S = SAFinit(M,mu,N,D,L);
tic;
[en,S] = SAFadapt(un,dn,S);

err_sqr = en.^2;

fprintf('Run time = %.2f s \n',toc);

figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:iter-1)/1024,10*log10(MSE), 'DisplayName', 'SAF');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'SAF' );

fprintf('\n');

%% MSAF
fprintf('MSAF\n');
fprintf('-------------------------------------------------------------\n');

S = MSAFinit(zeros(M,1),mu,N,L); % Initialization
tic;
[en,S] = MSAFadapt(un,dn,S);     % Perform MSAF algorithm

err_sqr = en.^2;

fprintf('Run time = %.2f s \n',toc);

figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:iter-1)/1024,10*log10(MSE), 'DisplayName', 'MSAF');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'MSAF' );

fprintf('\n');

%% SOAF
fprintf('SOAF\n');
fprintf('-------------------------------------------------------------\n');

S = SOAFinit(zeros(M,1),mu_small,iter);   % Initialization
tic;
[~,en,S] = SOAFadapt(un,dn,S);     % Perform algorithm

err_sqr = en.^2;

fprintf('Run time = %.2f s \n',toc);

figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:iter-1)/1024,10*log10(MSE), 'DisplayName', 'SOAF');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'SOAF' );

fprintf('\n');


%% Adding title and labels to plots
figure(MSEfig);
ylabel('Mean-square error'); grid on;
legend('show');
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)');


figure(NMSEfig);
axis([0 iter/1024 -20 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Cumulative Normalized Mean-square error'); grid on;
legend('show');


fprintf('\n');  % Empty line in logfile
