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

order = 1; 
M = 256; %% hammerstein filters lens


level = 2;                     % Levels of Wavelet decomposition
filters = 'db1';               % Set wavelet type
gains = rand(1,order)-0.5;
% gains = ones(1,order);

tube = @(x) tube(x,1); 
gdist = @(x) gdist(0.5,x); 
tanh = @(x) tanh(x)/2;  %from -0.5 to 0.5
pow = @(x) x.^3;
p = rand(1,order+1);
polynom = @(x) polyval(flip(p),x);


non_linearity = gdist; 
plot_nonlinearity = 1; 


iter = 0.5*80000;   % Number of iterations

%algorithm parameters
% p , w

mu_p = [0.5];
mu_w = [0.8];

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
NMSE_fig = figure('Name', 'NMSE');
nmse_n_points = 1000; 

%% RW data
filename = 'behr_ref_speech_m.wav'; %behr_ref_white  behr_ref_color  behr_ref_speech_f  behr_ref_speech_m
[dn, fs] = audioread(filename);
dn = dn(:,1)';

% import un variable
load opt_M.mat      %opt_white  opt_color  opt_M  opt_F
un = un';

% Resample the data to new fsn
fs_new = 8000;
[P, Q] = rat(fs_new/fs);
dn = resample(dn, P, Q);
un = resample(un, P, Q);

% % Normalization of u(n) and d(n)
% un = un/std(dn);
% dn = dn/std(dn);

% Normalization of speech signal
a = abs(max(dn))/sqrt(2);
un = un/a; dn = dn/a;

iter = length(un);

%% HAMMERSTERIN
fprintf('WAVPACK-HAMMERSTEIN\n');
fprintf('Order = %d, sys_len = %d \n', order, M);
for i = 1:runs
    fprintf('-------------------------------------------------------------\n');
    fprintf('Running iter (%d) of (%d)\n', i, runs);           
    fprintf('Run hyperpar: mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p(par_comb(1,i)), mu_w(par_comb(2,i)),alpha(par_comb(3,i)))


    [un,dn,vn] = GenerateResponses_nonlinear_Hammerstein(iter,lin_sys,sum(100*clock),4,40, non_linearity); %iter, b, seed, ARtype, SNR

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