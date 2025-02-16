% Volterra MWSAF       Multi Structured Wavelet-domain Adaptive Filter Demo
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

diary log_Volterra_TB.txt

fprintf('%s \n', datestr(datetime('now')));

addpath(genpath('../Common'));             % Functions in Common folder
addpath('Volterra NLMS fullband');
addpath('DFT_bank Volterra'); 
addpath('../MWSAF'); 
clear all;  
close all;

%% Hyperparameters
% Kernel Hyperpar
order = 2;                      % Order of volterra filter (just 2 atm)
M1 = 256;                       % length of first order volterra kernel
M2 = 32;                        % length of second order volterra kernel
M = [M1, M2];
gains = [1 1];                  % Kernel gains

% Signal Hyperpar
speech = 0;                     % Choose either 1 or 0, for using speech sample or noise 
% Choose ['speech_harvard_f.mat' ; 'speech_harvard_m.mat' ; 'SpeechSample.mat' ; 'speech.mat' ; 'Timit_m.mat' ; 'Timit_f_16k.mat' ] 
speech_sample = 'speech_harvard_m.mat';  
AR = 4;                         % AutoRegressive filter for white noise shaping, choose either 1 to 4, 1 is white noise
iter = 1*80000;                 % Number of iterations, NON speech
SNR = 40;

% WAVTERRA Hyperpar            (This can be modified to allow comparison)
par_level = [3];
par_filters = {'db16'};

par_C = M2;                     % Channels, #diagonal of kernel (max: M2)
par_SB = 1:2^par_level(end);    % Nonlinear subband (max: 1:2^level) NOT NEEDED

mu = [0.1, 0.1];                % Stepsize for different kernels 

% SAF & MSAF Hyperpar
N = 8;
D = N/2;                    % Decimation factor for 2x oversampling
L = 8*N;                    % Length of analysis filters, M=2KN, 

% Linear model hyperpar
l_mu = 0.1;
l_level = 3;
l_filters = 'db16';
Ml = M1;


% Create combination for Wavterra searchgrid
runs = length(par_level)*length(par_filters);
par_comb = combvec(1:length(par_level), 1:length(par_filters));

%% Create and plot kernel
% Create Kernel mode
kermode = 'simulated';              %modes: "delta", "randomdiag", "random", "lowpass", "simulated"

% Create Kernel parameters
deltapos = [5, 3];
Ndiag = 5;
normfreq = [0.6, 0.2];
h1h2 = ["h1.dat", "h2.dat"];
param = {deltapos, Ndiag, normfreq, h1h2};

[ker1, ker2] = create_kernel(M1, M2, kermode, param);

NL_system.M = [M1, M2];
NL_system.Responses = {gains(1).*ker1, gains(2).*ker2};

% Plotting kernel
kernel_plot(NL_system.Responses);

%% Create Signals
disp('Creating desired and input signals. . .');
if speech == 1
    fprintf('Kernel Memory: [%d, %d], Input signal: (%s) \n', M1, M2, speech_sample);
    [un,dn,vn] = GenerateResponses_speech_Volterra(NL_system, speech_sample);
    figure('Name', 'SpeechSpectrum');
    spectrogram(un, 1024, 256, 1024, 'yaxis');
else
    fprintf('Kernel Memory: [%d, %d], Input signal: noise colored wtih AR(%d)\n', M1, M2, AR);
    [un,dn,vn] = GenerateResponses_Volterra(iter, NL_system ,sum(100*clock),AR,SNR); %iter, b, seed, ARtype, SNR
    powerspec_plot(AR, iter);
end
fprintf('\n');


% figures handlers
MSEfig = figure('Name', 'MSE');
NMSEfig = figure('Name', 'NMSE');

%% WAVTERRA
fprintf('WAVTERRA\n');
for i = 1:runs  
    fprintf('--------------------------------------------------------------------\n'); 
    level = par_level(par_comb(1,i));
    filters = par_filters{par_comb(2,i)};
    SB = par_SB;
    C = par_C;
    
    fprintf('Running iter (%d) of (%d), level = %d , wtype = %s\n', i, runs, level, filters);           
    fprintf('step size = %s \n', sprintf('%.2f ', mu));

    tic;
    S = Volterra_Init(NL_system.M, mu, level, filters); 
    [en, S] = Volterra_2ord_adapt_v3(un, dn, S, C, SB);
 
    err_sqr = en.^2;
    er_len = length(err_sqr);

    fprintf('Total time = %.2f s \n',toc);

    % Plot MSE       
    figure(MSEfig);
    if speech == 1
        subplot(4, 1, 1:3);
    end
    q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
    hold on; plot((0:er_len-1)/1024,10*log10(MSE), 'DisplayName', ['WAVTERRA - Level:' ,num2str(level), ' ', filters]);
    
    NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
    fprintf('NMSE = %.2f dB\n', NMSE);
    
    % Plot NMSE
    figure(NMSEfig);
    hold on; plot((0:er_len-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['WAVTERRA - Level:' ,num2str(level), ' ', filters] );
            
    fprintf('\n');
end

%% MSAFTERRA
fprintf('MSAFTERRA\n');
fprintf('--------------------------------------------------------------------\n');
fprintf('Number of subbands, N = %d, step size = %s \n', N, sprintf('%.2f ', mu));

S = MSAFTERRA_Init(NL_system.M,mu,N,L);
tic;
[en,S] = MSAFTERRA_adapt(un,dn,S);
err_sqr = en.^2;

fprintf('Total time = %.2f s \n',toc);

figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:er_len-1)/1024,10*log10(MSE), 'DisplayName', ['MSAFTERRA, N: ', num2str(N), ' subband']);

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:er_len-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['MSAFTERRA, N: ', num2str(N), ' subband'] );

fprintf('\n');

%% SAFTERRA
fprintf('SAFTERRA\n');
fprintf('--------------------------------------------------------------------\n');
fprintf('Number of subbands, N = %d, Decimation factor, D = %d, step size = %s \n', N, D, sprintf('%.2f ', mu));

S = SAFTERRA_Init(NL_system.M,mu,N,D,L);
tic;
[en,S] = SAFTERRA_adapt(un,dn,S);
err_sqr = en.^2;

fprintf('Total time = %.2f s \n',toc);

figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:er_len-1)/1024,10*log10(MSE), 'DisplayName', ['SAFTERRA, N: ', num2str(N), ' subband']);

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:er_len-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['SAFTERRA, N: ', num2str(N), ' subband'] );

fprintf('\n');

%% FULLBAND VOLTERRA
fprintf('FULLBAND VOLTERRA NLMS\n');
fprintf('-------------------------------------------------------------\n');

Sfull = Volterra_NLMS_init(NL_system.M, mu); 

tic;
[en, Sfull] = Volterra_NLMS_adapt(un, dn, Sfull);     

err_sqr_full = en.^2;
    
fprintf('Total time = %.2f s \n',toc);

% Plot MSE
figure(MSEfig);
q = 0.99; MSE_full = filter((1-q),[1 -q],err_sqr_full);
plot((0:length(MSE_full)-1)/1024,10*log10(MSE_full), 'DisplayName', 'FB NLMS');

NMSE_FB = 10*log10(sum(err_sqr_full(256:end))/sum(dn(256:end).^2));
fprintf('NMSE = %.2f dB\n', NMSE_FB);

figure(NMSEfig);
hold on; plot((0:er_len-1)/1024, 10*log10(cumsum(err_sqr_full)./(cumsum(dn.^2))), 'DisplayName', 'FB NLMS');

fprintf('\n');


%% LINEAR MODEL
fprintf('LINEAR WMSAF\n');
fprintf('--------------------------------------------------------------------\n');
fprintf('Wavelet type: %s, levels: %d, step size = %.2f, filter length = %d\n', l_filters, l_level, l_mu, Ml);

tic;
Slin = SWAFinit(Ml, l_mu, l_level, l_filters); 
[en, Slin] = MWSAFadapt(un, dn, Slin); 

err_sqr = en.^2;

fprintf('Total time = %.2f s \n',toc);

figure(MSEfig);
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:er_len-1)/1024,10*log10(MSE), 'DisplayName', 'Linear WMSAF');

NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
fprintf('NMSE = %.2f dB\n', NMSE);

% Plot NMSE
figure(NMSEfig);
hold on; plot((0:er_len-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', 'Linear WMSAF' );

fprintf('\n');

%% Adding title and labels to plots
figure(MSEfig);
ylabel('Mean-square error'); grid on;
legend('show');
if speech == 1
    axis([0 er_len/1024 -120 inf]);
    subplot(4,1,4)
    plot((0:er_len-1)/1024, un, 'DisplayName', 'un');   
    axis([0 er_len/1024 -inf inf]);
    ylabel('Amplitude'); grid on;
else
    axis([0 er_len/1024 -60 10]);
end
xlabel('Number of iterations (\times 1024 input samples)');


figure(NMSEfig);
axis([0 er_len/1024 -20 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Cumulative Normalized Mean-square error'); grid on;
legend('show');


fprintf('\n');  % Empty line in logfile

diary off
