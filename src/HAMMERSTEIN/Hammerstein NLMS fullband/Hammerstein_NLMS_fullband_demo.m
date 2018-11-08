% Volterra NLMS fullband demo       Fullband Volterra second order adaptive
%                                   filtering
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

fprintf('%s \n', datestr(datetime('now')));
addpath(genpath('../../Common'));             % Functions in Common folder
clear all;  
close all;

%% Unidentified System parameters
<<<<<<< HEAD

order = 3; 
M = 32; %% hammerstein filters lens
gains = ones(3,1);

%algorithm parameters
mu = [0.02 0.02]; %ap aw
=======
order = 5; 
M = 256; %% hammerstein filters lens
gains = ones(5,1);

%algorithm parameters
mu = [0.1 0.5];
>>>>>>> 525178961d1e2ffd4a465bcb4e480049feeaf313
leak = [0 0]; 

%% Random Vector 
rng('default'); %

<<<<<<< HEAD


lin_system = load('h1.dat');
=======
lin_system = load('h2.dat');
>>>>>>> 525178961d1e2ffd4a465bcb4e480049feeaf313
lin_system = lin_system(1:M); 


%% Fullband Volterra NLMS
fprintf('-------------------------------------------------------------\n');
fprintf('FULLBAND VOLTERRA NLMS\n');

% Run parameters
iter = 4.0*80000;   % Number of iterations

<<<<<<< HEAD
[un,dn,vn] = GenerateResponses_nonlinear(iter,b,sum(100*clock),1,40, 'tanh', 0.2); %iter, b, seed, ARtype, SNR


%[un,dn,vn] = GenerateResponses_Hammerstein(iter, lin_system ,gains, order,sum(100*clock),1,40);
=======
[un,dn,vn] = GenerateResponses_Hammerstein(iter, lin_system ,gains, order,sum(100*clock),1,40);
>>>>>>> 525178961d1e2ffd4a465bcb4e480049feeaf313
% [un,dn,vn] = GenerateResponses_speech_Volterra(NL_system,'speech.mat');

tic;
Sfull = Hammerstein_NLMS_init(order, M, mu, leak); 
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
