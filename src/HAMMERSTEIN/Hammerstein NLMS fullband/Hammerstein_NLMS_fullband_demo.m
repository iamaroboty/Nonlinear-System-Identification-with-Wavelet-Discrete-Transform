% Volterra NLMS fullband demo       Fullband Volterra second order adaptive
%                                   filtering
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

fprintf('%s \n', datestr(datetime('now')));
addpath(genpath('../../Common'));             % Functions in Common folder
clear all;  
close all;

%% Unidentified System parameters
order = 3; 
M = 256; %% hammerstein filters lens
gains = ones(1,order);

iter = 0.5*80000;   % Number of iterations

%algorithm parameters
% p , w
leak = [0 0];
mu_p = [0.3 0.5 0.7];
mu_w = [0.3 0.5 0.7];
alpha = 10.^[0 -1 -2];

% Create combination
runs = length(mu_p)*length(mu_w)*length(alpha);
par_comb = combvec(1:length(mu_p), 1:length(mu_w),1:length(alpha));

% Random Vector 
rng('default'); %
lin_sys = rand(M,1);

% figures handlers
MSEfig = figure('Name', 'MSE');

%% HAMMERSTERIN
fprintf('HAMMERSTEIN\n');
fprintf('Order = %d, sys_len = %d \n', order, M);
for i = 1:runs
    fprintf('-------------------------------------------------------------\n');
    fprintf('Running iter (%d) of (%d)\n', i, runs);           
    fprintf('Run hyperpar: mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p(par_comb(1,i)), mu_w(par_comb(2,i)),alpha(par_comb(3,i)))

    [un,dn,vn] = GenerateResponses_Hammerstein(iter, lin_sys ,gains, order,sum(100*clock),1,40);

    tic;
    S = Hammerstein_NLMS_init(order, M, [mu_p(par_comb(1,i)) mu_w(par_comb(2,i))] ,leak, alpha(par_comb(3,i))); 

    [en, S] = Hammerstein_NLMS_adapt(un, dn, S);  

    err_sqr = en.^2;

    fprintf('Run time = %.2f s \n',toc);

    % Plot MSE
    figure(MSEfig);
    q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr);
    plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', ...
                                                            ['\mu_p:', num2str(mu_p(par_comb(1,i))),...
                                                             ' \mu_w:', num2str(mu_w(par_comb(2,i))),...
                                                             ' \alpha:', num2str(alpha(par_comb(3,i)))]);    
    axis([0 length(MSE_full)/1024 -inf 10]);
    xlabel('Number of iterations (\times 1024 input samples)'); 
    ylabel('Mean-square error (with delay)'); grid on; hold on;
    legend('show');

    NMSE(i) = 10*log10(sum(err_sqr)/sum(dn.^2));
    fprintf('NMSE = %.2f dB\n', NMSE(i));
end
    
[~, i] = min(NMSE);
fprintf('Hyperpar best sys: mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p(par_comb(1,i)), mu_w(par_comb(2,i)),alpha(par_comb(3,i)));

fprintf('\n');  % Empty line in logfile
