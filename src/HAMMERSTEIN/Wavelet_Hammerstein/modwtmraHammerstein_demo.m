% Volterra NLMS fullband demo       Fullband Volterra second order adaptive
%                                   filtering
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]


%% DFT Domain Hammerstein MSAF 

fprintf('%s \n', datestr(datetime('now')));
addpath(genpath('../../Common'));             % Functions in Common folder
clear all;  
close all;

%% Unidentified System parameters

M = 128; %% linear_part filt len

tube = @(x) tube(x,1); % tube distortion 
gdist = @(x) gdist(0.5,x); % guitar distotion 
tanh = @(x) tanh(x);  %from -0.5 to 0.5
pow = @(x) x.^2;
lin = @(x) x; 

% polinomial rep
order = 2; 
p = rand(1,order+1);
polynom = @(x) polyval(flip(p),x);

non_linearity = tube; 
plot_nonlinearity = 0; 

%% ALGORITHM PARAMETERS

mu = 0.01; 
palpha = 0.3;
alpha = 10^-4; 
wv = 'db4';

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
% input_type = {'white'};
% input_gain = {'75'};
% 
% fs = 44100;
% un = audioread([input_type{1},'_noise.wav']);
% % un = audioread('speech_f1_44100.wav');
% un = un(:)';
% dn = audioread(['univpm2_',input_type{1},'_',input_gain{1},'.wav']);
% %         dn = audioread(['univpm2_',input_type{i},'_',input_gain{j},'.wav']); 
% dn = dn(:)';
% 
% delay = 1024;
% dn = dn(delay:end); 
% un = un(1:end-delay+1); 
% 
% dn = dn(1:iter);
% un = un(1:iter);
% 
% % %     Resample the data to new fsn
% %         fs_new = 8000;
% %         [P, Q] = rat(fs_new/fs);
% %         dn = resample(dn, P, Q);
% %         un = resample(un, P, Q);
% %         iter = length(un);
% 
% %         Normalization of u(n) and d(n)
% %         un = un/std(dn);
% %         dn = dn/std(dn);
% un = un/(abs(max(un)));
% dn = dn/(abs(max(dn)));
% 
% 
% %         Normalization of speech signal
% %         a = abs(max(dn))/sqrt(2);
% %         un = un/a; dn = dn/a;

%% RUN WITH MULTIPLE PARAMETERS

% Create combination
runs = length(mu)*length(alpha);
par_comb = combvec(1:length(mu),1:length(alpha));

fprintf('WAVPACK-HAMMERSTEIN\n');
fprintf('Order = %d, sys_len = %d \n', order, M);
for i = 1:runs
    fprintf('-------------------------------------------------------------\n');
    fprintf('Running iter (%d) of (%d)\n', i, runs);           
    fprintf('Run hyperpar: mu = %.5f, alpha = %.2f \n', mu(par_comb(1,i)),alpha(par_comb(2,i)))
    
    
    %[un,dn,vn] = GenerateResponses_Hammerstein(iter, lin_sys , order, gains,sum(100*clock),2,40);
    %[un,dn,vn] = GenerateResponses_Hammerstein(iter, lin_sys , order, gains,sum(100*clock),2,40);    
    [un,dn,vn] = GenerateResponses_nonlinear_Hammerstein(iter,lin_sys,sum(100*clock),1,40, non_linearity); 

    tic; 
    
    S = modwtmraHammerstein_init(M, mu(par_comb(1,i)) , wv, palpha, alpha(par_comb(2,i))); 

    [en, S] = modwtmraHammerstein_adapt(un, dn, S);  

    err_sqr = en.^2;

    fprintf('Run time = %.2f s \n',toc);

    % Plot MSE
    figure(MSEfig);
    q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr);
    plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', ...
                                                            ['\mu:', num2str(mu(par_comb(1,i))),...                                                        
                                                             ' \alpha:', num2str(alpha(par_comb(2,i)))]);    
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
                    ['\mu:', num2str(mu(par_comb(1,i))),...
                     ' \alpha:', num2str(alpha(par_comb(2,i)))]);  
    axis([indx(1) indx(end) -15 5]);
    xlabel('Iteration number'); 
    ylabel('Cumulative Normalized Mean-square error'); 
    legend('show');   
end
    
[~, i] = min(NMSE);
fprintf('Hyperpar best sys: mu = %.5f, alpha = %.2f \n', mu(par_comb(1,i)), alpha(par_comb(2,i)));

fprintf('\n');  % Empty line in logfile

if plot_nonlinearity ==1 
    
    temp = real(ifft(S.weight));
    w = temp(1:M);
    
    figure;
    range = linspace(-10,10,1000); 
    plot(range,non_linearity(range)); 
    hold on; 
    poly = [0, S.coeffs{1,2}']; 
    ev = polyval(flip(poly), range); 
    plot(range, ev); 
    axis([min(range) max(range), -2, 2]); 
end