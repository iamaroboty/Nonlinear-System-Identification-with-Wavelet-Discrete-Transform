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

order = 5; 
M = 256; %% hammerstein filters lens
%gains = ones(8,1);

%algorithm parameters
mu = [0.02 0.02]; %ap aw
leak = [0 0]; 



un = load('ref_low_pass'); 
un = un.un;
dn = load('horn_low_pass'); 
dn = dn.dn;
 latency = 2880; 
 dn = dn(latency:end); 
 un = un(1:end-latency); 
max_iter = size(un,2); 



%% Fullband Volterra NLMS
fprintf('-------------------------------------------------------------\n');
fprintf('FULLBAND VOLTERRA NLMS\n');

% Run parameters
iter = 2.0*80000;   % Number of iterations

if iter > max_iter
   
    disp("WARNING: iter must be < max iter"); 
    iter = max_iter; 
    
end

un = un(1,1:iter); 
dn = dn(1,1:iter); 


%[un,dn,vn] = GenerateResponses_nonlinear(iter,b,sum(100*clock),1,40, 'tanh', 0.2); %iter, b, seed, ARtype, SNR


%[un,dn,vn] = GenerateResponses_Hammerstein(iter, lin_system ,gains, order,sum(100*clock),1,40);
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
diary off
