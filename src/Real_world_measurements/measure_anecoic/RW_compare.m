% Real World system identification from One Plus Two device
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]


addpath '../Common';             % Functions in Common folder

clear all;  
close all;
rng('default');

%% Structure Hyperparameters
% Filters length
M1 = 512;
M2 = 32;

% Universal stepsize
mu = [0.01 0.01];
mu_small = 0.003;

% For Wavelet Transform
level = 3;
filters = 'db8';

% Hammerstein and Wammerstein
% p , w
leak = [0 0];
order = 5;
M_hamm = M1 + M2;
mu_p = 0.3;
mu_w = 0.5;
alpha = 10.^0;

% For MSAF and SAF
N = 2^level;                % Number of subbands                     
D = N/2;                    % Decimation factor for 2x oversampling
L = 8*N;                    % Length of analysis filters, M=2KN

% Iter count
iter = 1.0*80000;

%input data 
infile = 'f'; %white, color, f, m
gain = '5p';  %5, 0, 5p


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
disp('Real-world desired and input signals. . .');
fprintf('Prior assumed Kernel Lengths: [%d, %d], iter= %d\n', M1, M2, iter);

% figures handlers
MSEfig = figure('Name', ['MSE (', infile, ' ', gain, ')']);
NMSEfig = figure('Name', ['NMSE (', infile, ' ', gain, ')']);

fprintf('\n');

%% WAVTERRA
fprintf('WAVTERRA\n');
fprintf('--------------------------------------------------------------------\n'); 
fprintf('level = %d , wtype = %s\n', level, filters);           
fprintf('step size = %s \n', sprintf('%.2f ', mu));

tic;
S = Volterra_Init([M1, M2], mu, level, filters); 
[en, S] = Volterra_2ord_adapt_v3(un, dn, S);

err_sqr = en.^2;

fprintf('Total time = %.2f s \n',toc);

% Plot MSE       
figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; grid on;
plot((1:iter)/1024,10*log10(MSE), 'DisplayName', ['WAVTERRA - Level:' ,num2str(level), ' ', filters]);
axis([0 iter/1024 -40 5]);
ylabel('Mean-square error'); grid on;
xlabel('Number of iterations (\times 1024 input samples)');
legend('show');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['WAVTERRA - Level:' ,num2str(level), ' ', filters] );
axis([0 iter/1024 -20 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Cumulative Normalized Mean-square error'); grid on;
legend('show');

fprintf('\n');

%% MSAFTERRA
% fprintf('MSAFTERRA\n');
% fprintf('--------------------------------------------------------------------\n');
% fprintf('Number of subbands, N = %d, step size = %s \n', N, sprintf('%.2f ', mu));
% 
% S = MSAFTERRA_Init([M1 M2],mu,N,L);
% tic;
% [en,S] = MSAFTERRA_adapt(un,dn,S);
% err_sqr = en.^2;
% 
% fprintf('Total time = %.2f s \n',toc);
% 
% % Plot MSE       
% figure(MSEfig);
% q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
% hold on; grid on;
% plot((1:iter)/1024,10*log10(MSE), 'DisplayName', ['MSAFTERRA, N: ', num2str(N), ' subband']);
% axis([0 iter/1024 -40 5]);
% ylabel('Mean-square error'); grid on;
% xlabel('Number of iterations (\times 1024 input samples)');
% legend('show');
% 
% NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
% fprintf('NMSE = %.2f dB\n', NMSE);
% 
% % Plot NMSE
% figure(NMSEfig);
% hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['MSAFTERRA, N: ', num2str(N), ' subband'] );
% axis([0 iter/1024 -20 10]);
% xlabel('Number of iterations (\times 1024 input samples)'); 
% ylabel('Cumulative Normalized Mean-square error'); grid on;
% legend('show');
% 
% fprintf('\n');


%% SAFTERRA
% fprintf('SAFTERRA\n');
% fprintf('--------------------------------------------------------------------\n');
% fprintf('Number of subbands, N = %d, Decimation factor, D = %d, step size = %s \n', N, D, sprintf('%.2f ', mu));
% 
% S = SAFTERRA_Init([M1 M2],mu,N,D,L);
% tic;
% [en,S] = SAFTERRA_adapt(un,dn,S);
% err_sqr = en.^2;
% 
% fprintf('Total time = %.2f s \n',toc);
% 
% figure(MSEfig);
% q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
% hold on; plot((1:iter)/1024,10*log10(MSE), 'DisplayName', ['SAFTERRA, N: ', num2str(N), ' subband'] );
% ylabel('Mean-square error'); grid on;
% legend('show');
% axis([0 iter/1024 -40 5]);
% xlabel('Number of iterations (\times 1024 input samples)');
% 
% NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
% fprintf('NMSE = %.2f dB\n', NMSE);
% 
% % Plot NMSE
% figure(NMSE_fig);  
% [NMSE, indx] = NMSE_compute(dn, en, nmse_n_points);
% hold on; grid on;
% plot(indx, NMSE,'DisplayName', ['SAFTERRA, N: ', num2str(N), ' subband']);
% axis([indx(1) indx(end) -10 0]);
% xlabel('Iteration number'); 
% ylabel('Cumulative Normalized Mean-square error'); 
% legend('show');
% 
% fprintf('\n');

%% FULLBAND VOLTERRA
fprintf('FULLBAND VOLTERRA NLMS\n');
fprintf('-------------------------------------------------------------\n');

S = Volterra_NLMS_init([M1 M2], mu); 

tic;
[en, S] = Volterra_NLMS_adapt(un, dn, S);     

err_sqr = en.^2;
    
fprintf('Total time = %.2f s \n',toc);

% Plot MSE       
figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; grid on;
plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'FB VOLTERRA');
axis([0 iter/1024 -40 5]);
ylabel('Mean-square error'); grid on;
xlabel('Number of iterations (\times 1024 input samples)');
legend('show');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'FB VOLTERRA' );
axis([0 iter/1024 -20 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Cumulative Normalized Mean-square error'); grid on;
legend('show');

fprintf('\n');


%% LINEAR MODEL
fprintf('LINEAR WMSAF\n');
fprintf('--------------------------------------------------------------------\n');
fprintf('Wavelet type: %s, levels: %d, step size = %.2f, filter length = %d\n', filters, level, mu(1), M1);

tic;
S = SWAFinit(M1, mu(1), level, filters); 
[en, S] = MWSAFadapt(un, dn, S ); 

err_sqr = en.^2;

fprintf('Total time = %.2f s \n',toc);

% Plot MSE       
figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; grid on;
plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'WMSAF');
axis([0 iter/1024 -40 5]);
ylabel('Mean-square error'); grid on;
xlabel('Number of iterations (\times 1024 input samples)');
legend('show');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'WMSAF' );
axis([0 iter/1024 -20 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Cumulative Normalized Mean-square error'); grid on;
legend('show');

fprintf('\n');

%% WAMMERSTERIN
fprintf('WAMMERSTERIN\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Order = %d, sys_len = %d \n', order, M_hamm);        
fprintf('mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p, mu_w,alpha)

tic;
S = Wammerstein_init(M_hamm, [mu_p mu_w] ,level, filters, order, alpha); 

[en, S] = Wammerstein_adapt(un, dn, S);  

err_sqr = en.^2;

fprintf('Run time = %.2f s \n',toc);

% Plot MSE       
figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; grid on;
plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'WAMMERSTEIN');
axis([0 iter/1024 -40 5]);
ylabel('Mean-square error'); grid on;
xlabel('Number of iterations (\times 1024 input samples)');
legend('show');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'WAMMERSTEIN' );
axis([0 iter/1024 -20 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Cumulative Normalized Mean-square error'); grid on;
legend('show');

fprintf('\n');

%% HAMMERSTEIN
fprintf('HAMMERSTEIN\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Order = %d, sys_len = %d \n', order, M_hamm);        
fprintf('mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p, mu_w,alpha)

tic;
S = Hammerstein_NLMS_init(order, M_hamm, [mu_p mu_w] ,leak, alpha); 

[en, S] = Hammerstein_NLMS_adapt(un, dn, S);  

err_sqr = en.^2;

fprintf('Run time = %.2f s \n',toc);

% Plot MSE       
figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; grid on;
plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'HAMMERSTEIN');
axis([0 iter/1024 -40 5]);
ylabel('Mean-square error'); grid on;
xlabel('Number of iterations (\times 1024 input samples)');
legend('show');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'HAMMERSTEIN' );
axis([0 iter/1024 -20 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Cumulative Normalized Mean-square error'); grid on;
legend('show');

fprintf('\n');


%% fix plots
if infile == 'm' || infile == 'f'
    figure(MSEfig);
    axis([0 iter/1024 -80 25]);
    
    figure(NMSEfig);
    axis([0 iter/1024 -40 20]);
end