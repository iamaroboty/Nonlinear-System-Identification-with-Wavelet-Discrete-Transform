% Volterra NLMS fullband demo       Fullband Volterra second order adaptive
%                                   filtering
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

diary log_VNLMS.txt
fprintf('%s \n', datestr(datetime('now')));
addpath(genpath('../../Common'));             % Functions in Common folder
clear all;  
close all;

%% Unidentified System parameters


NL_system.M = [256 16 2 4]; %% hammerstein filters lens
gains = [1 1 1 0];

%algorithm parameters
mu = [0.2, 0.2 0.01 0];
leaks = [0 0 0 0]; 

%% Random Vector 
rng('default'); %



NL_system = create_hammerstein_sys(NL_system.M, gains); 


%% Fullband Volterra NLMS
fprintf('-------------------------------------------------------------\n');
fprintf('FULLBAND VOLTERRA NLMS\n');

% Run parameters
iter = 3.0*80000;   % Number of iterations


[un,dn,vn] = GenerateResponses_Hammerstein(iter, NL_system ,sum(100*clock),1,40);
% [un,dn,vn] = GenerateResponses_speech_Volterra(NL_system,'speech.mat');

tic;
Sfull = Hammerstein_NLMS_init(NL_system.M, mu, leaks); 
[en, Sfull] = Hammerstein_NLMS_adapt(un, dn, Sfull);     

err_sqr_full = en.^2;
    
fprintf('Total time = %.3f mins \n',toc/60);

% Plot MSE
figure;
q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr_full);
plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', 'FB NLMS');
axis([0 length(MSE_full)/1024 -inf 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
legend('show');

fprintf('NMSE = %.2f dB\n', 10*log10(sum(err_sqr_full)/sum(dn.^2)));

fprintf('\n');  % Empty line in logfile
diary off
