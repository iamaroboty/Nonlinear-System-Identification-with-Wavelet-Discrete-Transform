% Real World system identification from One Plus Two device
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

diary opt_RW.txt

addpath '../Common';             % Functions in Common folder

clear all;  
close all;

%% Input Data
% filename = 'behr_ref_color.wav'; %behr_ref_white  behr_ref_color  behr_ref_speech_f  behr_ref_speech_m
% [dn, fs] = audioread(filename);
% dn = dn(:,1)';
% 
% % import un variable
% load opt_color.mat      %opt_white  opt_color  opt_M  opt_F
% un = un';

%Xperia values
load('xperia_ref_color.mat')
load('xperia_resp_color.mat')
delay = finddelay(un,dn);
dn = dn(abs(delay):end); 
un = un(1:end-abs(delay)+1); 

fs = 44100;

assert(isrow(un), 'un must be row vector');
assert(isrow(dn), 'dn must be row vector');
assert(length(dn) == length(un), 'ERROR in length of un and dn');

% Resample the data to new fsn
fs_new = 16000;
[P, Q] = rat(fs_new/fs);
dn = resample(dn, P, Q);
un = resample(un, P, Q);

% Normalization of u(n) and d(n)
un = un/std(dn);
dn = dn/std(dn);

% % Normalization of speech signal
% a = abs(max(dn))/sqrt(2);
% un = un/a; dn = dn/a;

noise_level = -120; % dBfs [measured] not real

%% Structure Hyperparameters
% Filters length
M1 = 512;
M2 = 32;

% For Wavelet Transform
level = 2;
filters = 'db8';

% For Volterra Based Methods
mu = [0.3 0.3];
C = M2;                     % Number of diagonals
SB = 1:2^level;             % Nonlinear Subbands

% For MSAFTERRA and SAFTERRA
N = 2^level;                % Number of subbands                     
D = N/2;                    % Decimation factor for 2x oversampling
L = 8*N;                    % Length of analysis filters, M=2KN

% Hammerstein Full Band
% p , w
order = 3;
M = M1;
leak = [0 0];
mu_p_fb = 0.3;
mu_w_fb = 0.5;
alpha_fb = 10.^2;

%Wammerstein
mu_p = 0.3;
mu_w = 0.5;
alpha = 10.^2;


% Iter count
iter = length(dn);

%%
disp('Real-world desired and input signals. . .');
fprintf('Prior assumed Kernel Lengths: [%d, %d], iter= %d\n', M1, M2, iter);

% figures handlers
MSEfig = figure('Name', 'MSE');
NMSE_fig = figure('Name', 'NMSE');
nmse_n_points = 1000; 

fprintf('\n');

%% WAVTERRA
% fprintf('WAVTERRA\n');
% fprintf('--------------------------------------------------------------------\n'); 
% fprintf('level = %d , wtype = %s\n', level, filters);           
% fprintf('step size = %s \n', sprintf('%.2f ', mu));
% 
% tic;
% S = Volterra_Init([M1, M2], mu, level, filters); 
% [en, S] = Volterra_2ord_adapt_v3(un, dn, S);
% 
% err_sqr = en.^2;
% 
% fprintf('Total time = %.2f s \n',toc);
% 
% % Plot MSE       
% figure(MSEfig);
% q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
% hold on; grid on;
% plot((1:iter)/1024,10*log10(MSE), 'DisplayName', ['WAVTERRA - Level:' ,num2str(level), ' ', filters]);
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
% plot(indx, NMSE,'DisplayName', ['WAVTERRA - Level:' ,num2str(level), ' ', filters]);
% axis([indx(1) indx(end) -10 0]);
% xlabel('Iteration number'); 
% ylabel('Cumulative Normalized Mean-square error'); 
% legend('show');
% 
% fprintf('\n');

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
% figure(MSEfig);
% q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
% hold on; plot((1:iter)/1024,10*log10(MSE), 'DisplayName', ['MSAFTERRA, N: ', num2str(N), ' subband']);
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
% plot(indx, NMSE,'DisplayName',  ['MSAFTERRA, N: ', num2str(N), ' subband']);
% axis([indx(1) indx(end) -10 0]);
% xlabel('Iteration number'); 
% ylabel('Cumulative Normalized Mean-square error'); 
% legend('show');
% 
% fprintf('\n');
% 
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
% fprintf('FULLBAND VOLTERRA NLMS\n');
% fprintf('-------------------------------------------------------------\n');
% 
% S = Volterra_NLMS_init([M1 M2], mu); 
% 
% tic;
% [en, S] = Volterra_NLMS_adapt(un, dn, S);     
% 
% err_sqr = en.^2;
%     
% fprintf('Total time = %.2f s \n',toc);
% 
% % Plot MSE
% figure(MSEfig);
% q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
% plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'FB NLMS');
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
% plot(indx, NMSE,'DisplayName', 'FB NLMS');
% axis([indx(1) indx(end) -10 0]);
% xlabel('Iteration number'); 
% ylabel('Cumulative Normalized Mean-square error'); 
% legend('show');
% 
% fprintf('\n');


%% LINEAR MODEL
fprintf('LINEAR WMSAF\n');
fprintf('--------------------------------------------------------------------\n');
fprintf('Wavelet type: %s, levels: %d, step size = %.2f, filter length = %d\n', filters, level, mu(1), M1);

tic;
S = SWAFinit(M1, mu(1), level, filters); 
[en, S ] = MWSAFadapt(un, dn, S ); 

err_sqr = en.^2;

fprintf('Total time = %.2f s \n',toc);

figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'Linear WMSAF');
ylabel('Mean-square error'); grid on;
legend('show');
axis([0 iter/1024 -40 5]);
xlabel('Number of iterations (\times 1024 input samples)');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSE_fig);  
[NMSE, indx] = NMSE_compute(dn, en, nmse_n_points);
hold on; grid on;
plot(indx, NMSE,'DisplayName', 'Linear WMSAF');
axis([indx(1) indx(end) -10 0]);
xlabel('Iteration number'); 
ylabel('Cumulative Normalized Mean-square error'); 
legend('show');

fprintf('\n');

%% WAMMERSTERIN
fprintf('WAMMERSTERIN\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Order = %d, sys_len = %d \n', order, M);        
fprintf('mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p, mu_w,alpha)

tic;
S = Wammerstein_init(M, [mu_p mu_w] ,level, filters, order, alpha); 

[en, S] = Wammerstein_adapt(un, dn, S);  

err_sqr = en.^2;

fprintf('Run time = %.2f s \n',toc);

% Plot MSE
figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr); hold on;
plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'WAMMERSTEIN');    
ylabel('Mean-square error'); grid on;
legend('show');
axis([0 iter/1024 -40 5]);
xlabel('Number of iterations (\times 1024 input samples)');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSE_fig);  
[NMSE, indx] = NMSE_compute(dn, en, nmse_n_points);
hold on; grid on;
plot(indx, NMSE,'DisplayName',  'WAMMERSTEIN');
axis([indx(1) indx(end) -10 0]);
xlabel('Iteration number'); 
ylabel('Cumulative Normalized Mean-square error'); 
legend('show');

fprintf('\n');

%% HAMMERSTEIN
fprintf('HAMMERSTEIN\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Order = %d, sys_len = %d \n', order, M);        
fprintf('mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p_fb, mu_w_fb,alpha_fb)

tic;
S = Hammerstein_NLMS_init(order, M, [mu_p mu_w] ,leak, alpha); 

[en, S] = Hammerstein_NLMS_adapt(un, dn, S);  

err_sqr = en.^2;

fprintf('Run time = %.2f s \n',toc);

% Plot MSE
figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr); hold on;
plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'HAMMERSTEIN');    
ylabel('Mean-square error'); grid on;
legend('show');
axis([0 iter/1024 -40 5]);
xlabel('Number of iterations (\times 1024 input samples)');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSE_fig);  
[NMSE, indx] = NMSE_compute(dn, en, nmse_n_points);
hold on; grid on;
plot(indx, NMSE,'DisplayName','HAMMERSTEIN');
axis([indx(1) indx(end) -10 0]);
xlabel('Iteration number'); 
ylabel('Cumulative Normalized Mean-square error'); 
legend('show');

fprintf('\n');

diary off
