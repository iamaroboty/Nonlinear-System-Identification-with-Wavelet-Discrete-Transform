% Volterra NLMS fullband demo       Fullband Volterra second order adaptive
%                                   filtering
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]


%% DFT Domain Hammerstein MSAF 

fprintf('%s \n', datestr(datetime('now')));
addpath(genpath('../../Common'));             % Functions in Common folder
clear all;  
close all;

%% run parameters 
rng('default'); % fix random geneerator 


%% Unidentified System parameters

M = 256; %% linear_part filt len




tube = @(x) tube(x,1); % tube distortion 
gdist = @(x) gdist(0.5,x); % guitar distotion 
tanh = @(x) tanh(x)/2;  %from -0.5 to 0.5
pow = @(x) x.^2;
lin = @(x) x; 

% polinomial rep
order = 1; 
p = rand(1,order+1);
polynom = @(x) polyval(flip(p),x);


non_linearity = tanh; 
plot_nonlinearity = 1; 

%% ALGORITHM PARAMETERS

order = 4; 

M_alg = 64; 

N = 4;          % Number of Subbands
D = N/2;        % Decimation factor for 2x oversampling
L = 8*N;        % Length of analysis filters, M=2KN, 


%learning rates for 
% p , w

mu_p = [0.2]; % taylor part
mu_w = [0.5]; % linear part

alpha = 10.^[1];

iter = 0.5*80000;   % Number of iterations


%% LINEAR PART SYSTEM PARAMETERS

% Random Vector 
b = load("h1.dat"); % load real world IR
b = b(1:M);
% lin_sys = rand(M,1)-0.5; % random IR
lin_sys = b;


% figures handlers
MSEfig = figure('Name', 'MSE');
NMSE_fig = figure('Name', 'NMSE');
nmse_n_points = 1000; 

%% RW data (if used set here and comment in the next section generatehammerstein fucntions and similar)

% un = audioread('ar10_noise.wav')';
% dn = audioread('univpm2_ar10_75.wav')'; 
% dn = dn(1:size(un,2)); 
% latency = 1024; 
% un = un(1:end-latency); 
% dn = dn(latency:end); 

% Resample the data to new fsn

% fs_new = 8000;
% [P, Q] = rat(fs_new/44100);
% dn = resample(dn, P, Q);
% un = resample(un, P, Q);

% % Normalization of u(n) and d(n)
% un = un/std(dn);
% dn = dn/std(dn);

% Normalization of speech signal
% a = abs(max(dn))/sqrt(2);
% un = un/a; dn = dn/a;
% 
% max_iter = length(un); 
% 
% if iter > max_iter
%    
%    iter = max_iter;     
% end

%% RUN WITH MULTIPLE PARAMETERS

% Create combination
runs = length(mu_p)*length(mu_w)*length(alpha);
par_comb = combvec(1:length(mu_p), 1:length(mu_w),1:length(alpha));

fprintf('WAVPACK-HAMMERSTEIN\n');
fprintf('Order = %d, sys_len = %d \n', order, M);
for i = 1:runs
    fprintf('-------------------------------------------------------------\n');
    fprintf('Running iter (%d) of (%d)\n', i, runs);           
    fprintf('Run hyperpar: mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p(par_comb(1,i)), mu_w(par_comb(2,i)),alpha(par_comb(3,i)))
    
    
    %[un,dn,vn] = GenerateResponses_Hammerstein(iter, lin_sys , order, gains,sum(100*clock),2,40);
    %[un,dn,vn] = GenerateResponses_Hammerstein(iter, lin_sys , order, gains,sum(100*clock),2,40);    
    [un,dn,vn] = GenerateResponses_nonlinear_Hammerstein(iter,lin_sys,sum(100*clock),1,40, non_linearity); 

    tic; 
    
    S = DFThammerstein_init(M, [mu_p(par_comb(1,i)) mu_w(par_comb(2,i))] , order, alpha(par_comb(3,i)), N, D , L); 

    [en, S] = DFThammerstein_adapt(un, dn, S);  

    err_sqr = en.^2;

    fprintf('Run time = %.2f s \n',toc);

    % Plot MSE
    figure(MSEfig);
    q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr);
    plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', ...
                                                            ['\mu_p:', num2str(mu_p(par_comb(1,i))),...
                                                             ' \mu_w:', num2str(mu_w(par_comb(2,i))),...
                                                             ' \alpha:', num2str(alpha(par_comb(3,i)))]);    
    axis([0 length(MSE_full)/1024 -50 10]);
    xlabel('Number of iterations (\times 1024 input samples)'); 
    ylabel('Mean-square error (with delay)'); grid on; hold on;
    legend('show');

    NMSE(i) = 10*log10(sum(err_sqr)/sum(dn.^2));
    fprintf('NMSE = %.2f dB\n', NMSE(i));
    
    figure(NMSE_fig);  
    [NMSE_plot, indx] = NMSE_compute(dn, en, nmse_n_points);
    hold on; grid on;
    plot(indx, NMSE_plot,      	  'DisplayName', ...
                    ['\mu_p:', num2str(mu_p(par_comb(1,i))),...
                     ' \mu_w:', num2str(mu_w(par_comb(2,i))),...
                     ' \alpha:', num2str(alpha(par_comb(3,i)))]);  
    axis([indx(1) indx(end) -15 5]);
    xlabel('Iteration number'); 
    ylabel('Cumulative Normalized Mean-square error'); 
    legend('show');   
end
    
[~, i] = min(NMSE);
fprintf('Hyperpar best sys: mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p(par_comb(1,i)), mu_w(par_comb(2,i)),alpha(par_comb(3,i)));

fprintf('\n');  % Empty line in logfile

if plot_nonlinearity ==1 
    figure;
    range = linspace(-10,10,1000); 
    plot(range,non_linearity(range)); 
    hold on; 
    poly = [0, S.coeffs{1,2}']; 
    ev = polyval(flip(poly), range); 
    plot(range, ev); 
    axis([min(range) max(range), -2, 2]); 
end