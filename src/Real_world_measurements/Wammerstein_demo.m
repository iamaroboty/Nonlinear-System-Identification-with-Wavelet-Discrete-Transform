% Volterra NLMS fullband demo       Fullband Volterra second order adaptive
%                                   filtering
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

fprintf('%s \n', datestr(datetime('now')));
addpath(genpath('../../Common'));             % Functions in Common folder
clear all;  
% close all;
rng('default'); %

%% Unidentified System parameters
order = 10; 
M = 256; %% hammerstein filters lens
level = 3;                     % Levels of Wavelet decomposition
filters = 'db4';               % Set wavelet type
gains = rand(1,order)-0.5;
% gains = ones(1,order);


iter = 2*80000;   % Number of iterations

%% dn signal

un = load('guitar_in'); 
un = un.un;
dn = load('engl_5'); 
dn = dn.dn;

[P,Q] = rat(8192/44100);

un = resample(un,P,Q); 
dn = resample(dn,P,Q); 

%estimate latency with xcorr
[corr, lag]= xcorr(un,dn); 
[~, ind]= max(abs(corr)); 
latency = abs(lag(ind)); 

% latency = 3130; %horn 2800 behr 3130 xperia 
dn = dn(latency:end); 
un = un(1:end-latency); 

% normalization 

% un = un/std(dn);
% dn = dn/std(dn);

%normalization speech 
a = abs(max(dn))/sqrt(2); 
un = un/a; 
dn= dn/ a; 

max_iter = size(un,2); 

if iter > max_iter
   
    disp("WARNING: iter must be < max iter"); 
    iter = max_iter; 
    
end

un = un(1,1:iter); 
dn = dn(1,1:iter); 


%algorithm parameters
% p , w
mu_p = [0.3];
mu_w = [0.7];
alpha = 10.^[1];


% Create combination
runs = length(mu_p)*length(mu_w)*length(alpha);
par_comb = combvec(1:length(mu_p), 1:length(mu_w),1:length(alpha));

% Random Vector 
b = load("h1.dat");
b = b(1:M);
% lin_sys = rand(M,1);
lin_sys = b;

% non_linearity = 'tanh'; %pow %tanh %relu %logsig
% fprintf('Non Linearity = %s \n', non_linearity);

% figures handlers
MSEfig = figure('Name', 'MSE');

%% HAMMERSTERIN
fprintf('WAVPACK-HAMMERSTEIN\n');
fprintf('Order = %d, sys_len = %d \n', order, M);
for i = 1:runs
    fprintf('-------------------------------------------------------------\n');
    fprintf('Running iter (%d) of (%d)\n', i, runs);           
    fprintf('Run hyperpar: mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p(par_comb(1,i)), mu_w(par_comb(2,i)),alpha(par_comb(3,i)))

     %[un,dn,vn] = GenerateResponses_Hammerstein(iter, lin_sys , order, gains,sum(100*clock),2,40);
    
    %[un,dn,vn] = GenerateResponses_nonlinear_Hammerstein(iter,lin_sys,sum(100*clock),4,40, non_linearity); %iter, b, seed, ARtype, SNR

    tic;
    S = Wammerstein_init(M, [mu_p(par_comb(1,i)) mu_w(par_comb(2,i))] ,level, filters, order, alpha(par_comb(3,i))); 

    [en, S] = Wammerstein_adapt(un, dn, S);  

    err_sqr = en.^2;

    fprintf('Run time = %.2f s \n',toc);

    % Plot MSE
    figure(MSEfig);
    q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr);
    subplot(2,1,1);
    plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', ...
                                                            ['\mu_p:', num2str(mu_p(par_comb(1,i))),...
                                                             ' \mu_w:', num2str(mu_w(par_comb(2,i))),...
                                                             ' \alpha:', num2str(alpha(par_comb(3,i)))]);    
    axis([0 length(MSE_full)/1024 -150 10]);
    xlabel('Number of iterations (\times 1024 input samples)'); 
    ylabel('Mean-square error (with delay)'); grid on; hold on;
    legend('show');

    NMSE(i) = 10*log10(sum(err_sqr)/sum(dn.^2));
    fprintf('NMSE = %.2f dB\n', NMSE(i));
end
    
[~, i] = min(NMSE);
fprintf('Hyperpar best sys: mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p(par_comb(1,i)), mu_w(par_comb(2,i)),alpha(par_comb(3,i)));

fprintf('\n');  % Empty line in logfile

subplot(2,1,2); 
plot((0:length(dn)-1)/1024,dn, 'DisplayName', 'dn');
axis([0 length(dn)/1024 -3  3]);
