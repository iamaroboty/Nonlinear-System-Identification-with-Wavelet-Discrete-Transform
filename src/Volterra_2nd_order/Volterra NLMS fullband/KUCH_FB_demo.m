% Volterra NLMS fullband demo       Fullband Volterra second order adaptive
%                                   filtering
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

diary log_Kuch.txt
fprintf('%s \n', datestr(datetime('now')));
addpath(genpath('../../Common'));             % Functions in Common folder
clear all;  
close all;

%% Unidentified System parameters
order = 2; 
M1 = 256; % length of first order volterra kernel
M2 = 32; % length of second order volterra kernel

NL_system.M = [M1, M2];
gains = [1 1];

%NL_system = create_volterra_sys(order, M, gains, 'nlsys1'); 


%% Random Vector 
rng('default'); % For reproducibility
kermode = 'randomdiag';     %modes: "delta", "randomdiag", "random", "lowpass", "simulated"

% Create Kernel parameters
deltapos = [1, 1];
Ndiag = 5;
normfreq = [0.6, 0.2];
h1h2 = ["h1.dat", "h2.dat"];
param = {deltapos, Ndiag, normfreq, h1h2};

[ker1, ker2] = create_kernel(M1, M2, kermode, param);

NL_system.M = [M1, M2];
NL_system.Responses = {gains(1).*ker1, gains(2).*ker2};

kernel_plot(NL_system.Responses);

% Run parameters
iter = 1.0*80000;                % Number of iterations

[un,dn,vn] = GenerateResponses_Volterra(iter, NL_system ,sum(100*clock),1,40);
% [un,dn,vn] = GenerateResponses_speech_Volterra(NL_system,'speech.mat');

%% Kuch
fprintf('-------------------------------------------------------------\n');
fprintf('Kuch\n');

mu = [0.1, 0.001, 0.001];       % mu1 muc muw
Nc = 15;

posdelta = 2;
Sfull = KUCH_Init(NL_system.M, Nc, mu, posdelta); 
tic;
[en, Sfull] = KUCH_NLMS_adapt(un, dn, Sfull);     

err_sqr_full = en.^2;
    
fprintf('Total time = %.2f s \n',toc);

% Plot MSE
figure;
q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr_full); hold on;
plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', ['Kuch, Nc = ', num2str(Nc)]);
axis([0 length(MSE_full)/1024 -inf 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
legend('show');

fprintf('NMSE = %.2f dB\n', 10*log10(sum(err_sqr_full)/sum(dn.^2)));
fprintf('\n');

% 
% %% Fullband Volterra NLMS
% fprintf('-------------------------------------------------------------\n');
% fprintf('FULLBAND VOLTERRA NLMS\n');
% 
% mu = [0.1, 0.1];
% 
% Sfull = Volterra_NLMS_init(NL_system.M, mu); 
% tic;
% [en, Sfull] = Volterra_NLMS_adapt(un, dn, Sfull);     
% 
% err_sqr_full = en.^2;
%     
% fprintf('Total time = %.2f s \n',toc);
% 
% % Plot MSE
% q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr_full);
% plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', 'FB NLMS');
% axis([0 length(MSE_full)/1024 -inf 10]);
% xlabel('Number of iterations (\times 1024 input samples)'); 
% ylabel('Mean-square error (with delay)'); grid on;
% legend('show');
% 
% fprintf('NMSE = %.2f dB\n', 10*log10(sum(err_sqr_full)/sum(dn.^2)));
% 
% 
fprintf('\n');  % Empty line in logfile
diary off



