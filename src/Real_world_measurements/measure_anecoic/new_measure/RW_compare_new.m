% Real World system identification from One Plus Two device
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]


addpath '../Common';             % Functions in Common folder

clear all;  
close all;

%% Structure Hyperparameters
% Filters length
M1 = 256;
M2 = 128;

% Universal stepsize
mu = [0.1 0.3];
mu_fb = [0.1 0.1];

% For Wavelet Transform
level = 3;
filters = 'db16';

% For Waveterra
C = 10;
SB = 1:2^level;

% Hammerstein and Wammerstein
% p , w
leak = [0 0];
order = 5;
M_hamm = M1 + M2;
mu_p = 0.3;
mu_w = 0.7;
alpha = 10.^0;

% For MSAF and SAF
N = 2^level;                % Number of subbands                     
D = N/2;                    % Decimation factor for 2x oversampling
L = 8*N;                    % Length of analysis filters, M=2KN

% Iter count
iter = 2*80000;

input_type = {'white'};
input_gain = {'75'};

n = 0;
for i= 1:numel(input_type)
    for j= 1:numel(input_gain)
        n = n+1;
        fprintf('STARTING RUN (%d)/(%d), file: %s \n', n , numel(input_type)*numel(input_gain), [input_type{i}, '_', input_gain{j}]);
        %% Input Data
        fs = 44100;
%         un = audioread([input_type{i},'_noise.wav']);
        un = audioread('speech_f1_44100.wav');
        un = un(:)';
        dn = audioread(['univpm2_',input_type{i},'_',input_gain{j},'.wav']);
%         dn = audioread(['univpm2_',input_type{i},'_',input_gain{j},'.wav']); 
        dn = dn(:)';
        
        delay = 1024;
        dn = dn(delay:end); 
        un = un(1:end-delay+1); 

        dn = dn(1:iter);
        un = un(1:iter);
 
% %     Resample the data to new fsn
%         fs_new = 8000;
%         [P, Q] = rat(fs_new/fs);
%         dn = resample(dn, P, Q);
%         un = resample(un, P, Q);
%         iter = length(un);
        
%         Normalization of u(n) and d(n)
%         un = un/std(dn);
%         dn = dn/std(dn);
        un = un/(abs(max(un)));
        dn = dn/(abs(max(dn)));
        

%         Normalization of speech signal
%         a = abs(max(dn))/sqrt(2);
%         un = un/a; dn = dn/a;

        %% Create and plot kernel simulated
%         % Create Kernel mode
%         kermode = 'simulated';              %modes: "delta", "randomdiag", "random", "lowpass", "simulated"
% 
%         % Create Kernel parameters
%         deltapos = [5, 3];
%         Ndiag = 1;
%         normfreq = [0.6, 0.2];
%         h1h2 = ["h1.dat", "h2.dat"];
%         param = {deltapos, Ndiag, normfreq, h1h2};
% 
%         [ker1, ker2] = create_kernel(M1, M2, kermode, param);
% 
%         NL_system.M = [M1, M2];
%         gains = [1 0.1]; 
%         NL_system.Responses = {gains(1).*ker1, gains(2).*ker2};
% 
%         % Plotting kernel
%         figure('Name', 'Simulated Kernel');
%         kernel_plot(NL_system.Responses);
%         [un,dn,vn] = GenerateResponses_Volterra(iter, NL_system ,sum(100*clock),4,40); %iter, b, seed, ARtype, SNR
%         iter = length(un);
%         % b = load('h1.dat');
%         % b = b(1:M);
%         % [un,dn,vn] = GenerateResponses(iter,b,sum(100*clock),1,40); %iter, b, seed, ARtype, SNR
%         % [un,dn,vn] = GenerateResponses_speech(b,'SpeechSample.mat');
        
        %%
        disp('Real-world desired and input signals. . .');
        fprintf('Prior assumed Kernel Lengths: [%d, %d], iter= %d\n', M1, M2, iter);

        % figures handlers
        % MSEfig = figure('Name', ['MSE (', infile, ' ', gain, ')']);
        % NMSEfig = figure('Name', ['NMSE (', infile, ' ', gain, ')']);
        fig = figure('Name', [input_type{i}, '_', input_gain{j}]);

        fprintf('\n');

        %% WAVTERRA
        fprintf('WAVTERRA\n');
        fprintf('--------------------------------------------------------------------\n'); 
        fprintf('level = %d , wtype = %s\n', level, filters);           
        fprintf('diags = %d, SB = %d, step size = %s \n', C, numel(SB), sprintf('%.2f ', mu));
        
        tic;
        S = Volterra_Init([M1, M2], mu, level, filters); 
        [en, S] = Volterra_2ord_adapt_v3(un, dn, S, C, SB);
        
        coeffs.wavterra = S.coeffs;
        figure('Name', ['Estimated kernel Wavterra ', input_type{i}, ' ', input_gain{j}]);
        visualize_est_ker(coeffs.wavterra);
%         savefig(['figures/Ker_WAVTERRA_',input_type{i},'_',input_gain{j},'.fig'])

        
        err_sqr = en.^2;
        
        fprintf('Total time = %.2f s \n',toc);
        
        % Plot MSE       
        figure(fig);
        subplot(211);
        q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
        hold on; grid on;
        plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'WAVTERRA');
        axis([0 iter/1024 -40 5]);
        legend('show');
        
        NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
        fprintf('NMSE = %.2f dB\n', NMSE);
        
        % Plot NMSE
        subplot(212);
        hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['WAVTERRA: ', num2str(num2str(NMSE, '%0.2f'))] );
        axis([0 iter/1024 -20 10]);
        legend('show');
        
        fprintf('\n');

        %% MSAFTERRA
        fprintf('MSAFTERRA\n');
        fprintf('--------------------------------------------------------------------\n');
        fprintf('Number of subbands, N = %d, step size = %s \n', N, sprintf('%.2f ', mu));
        
        S = MSAFTERRA_Init([M1 M2],mu,N,L);
        tic;
        [en,S] = MSAFTERRA_adapt(un,dn,S, C);
        coeffs.msafterra = S.coeffs;
        err_sqr = en.^2;
        
        fprintf('Total time = %.2f s \n',toc);
        
        % Plot MSE       
        figure(fig);
        subplot(211);
        q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
        hold on; grid on;
        plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'MSAFTERRA');
        axis([0 iter/1024 -40 5]);
        legend('show');
        
        NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
        fprintf('NMSE = %.2f dB\n', NMSE);
        
        % Plot NMSE
        subplot(212);
        hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['MSAFTERRA: ', num2str(num2str(NMSE, '%0.2f'))] );
        axis([0 iter/1024 -20 10]);
        legend('show');
        
        fprintf('\n');

        %% SAFTERRA
%         fprintf('SAFTERRA\n');
%         fprintf('--------------------------------------------------------------------\n');
%         fprintf('Number of subbands, N = %d, Decimation factor, D = %d, step size = %s \n', N, D, sprintf('%.2f ', mu));
%         
%         S = SAFTERRA_Init([M1 M2],mu,N,D,L);
%         tic;
%         [en,S] = SAFTERRA_adapt(un,dn,S, C);
%         coeffs.safterra = S.coeffs;
%         err_sqr = en.^2;
%         
%         fprintf('Total time = %.2f s \n',toc);
%         
%         % Plot MSE       
%         figure(fig);
%         subplot(211);
%         q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
%         hold on; grid on;
%         plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'SAFTERRA');
%         axis([0 iter/1024 -40 5]);
%         legend('show');
%         
%         NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
%         fprintf('NMSE = %.2f dB\n', NMSE);
%         
%         % Plot NMSE
%         subplot(212);
%         hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['SAFTERRA: ', num2str(num2str(NMSE, '%0.2f'))] );
%         axis([0 iter/1024 -20 10]);
%         legend('show');
%         
%         fprintf('\n');

        %% FULLBAND VOLTERRA
        fprintf('FULLBAND VOLTERRA NLMS\n');
        fprintf('-------------------------------------------------------------\n');

        S = Volterra_NLMS_init([M1 M2], mu_fb); 

        tic;
        [en, S] = Volterra_NLMS_adapt(un, dn, S);  
        
        coeffs.volterra = S.coeffs;
        figure('Name', ['Estimated kernel Volterra ', input_type{i}, ' ', input_gain{j}]);
        visualize_est_ker(coeffs.volterra);
%         savefig(['figures/Ker_Volterra_',input_type{i},'_',input_gain{j},'.fig'])
        
        err_sqr = en.^2;

        fprintf('Total time = %.2f s \n',toc);

        % Plot MSE       
        figure(fig);
        subplot(211);
        q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
        hold on; grid on;
        plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'FB VOLTERRA');
        axis([0 iter/1024 -40 5]);
        legend('show');

        NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
        fprintf('NMSE = %.2f dB\n', NMSE);

        % Plot NMSE
        subplot(212);
        hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['FB VOLTERRA: ', num2str(num2str(NMSE, '%0.2f'))] );
        axis([0 iter/1024 -20 10]);
        legend('show');

        fprintf('\n');


        %% LINEAR MODEL
        fprintf('LINEAR WMSAF\n');
        fprintf('--------------------------------------------------------------------\n');
        fprintf('Wavelet type: %s, levels: %d, step size = %.2f, filter length = %d\n', filters, level, mu(1), M1);

        tic;
        S = SWAFinit(M1, mu(1), level, filters); 
        [en, S] = MWSAFadapt(un, dn, S ); 
        coeffs.wmsaf = S.coeffs;

        err_sqr = en.^2;

        fprintf('Total time = %.2f s \n',toc);

        % Plot MSE       
        figure(fig);
        subplot(211);
        q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
        hold on; grid on;
        plot((1:iter)/1024,10*log10(MSE), 'DisplayName', 'WMSAF');
        axis([0 iter/1024 -40 5]);
        legend('show');

        NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
        fprintf('NMSE = %.2f dB\n', NMSE);

        % Plot NMSE
        subplot(212);
        hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['WMSAF: ', num2str(num2str(NMSE, '%0.2f'))] );
        axis([0 iter/1024 -20 10]);
        legend('show');

        fprintf('\n');

        %% WAMMERSTERIN
        fprintf('WAMMERSTERIN\n');
        fprintf('-------------------------------------------------------------\n');
        fprintf('Order = %d, sys_len = %d \n', order, M_hamm);        
        fprintf('mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p, mu_w,alpha)

        tic;
        S = Wammerstein_init(M_hamm, [mu_p mu_w] ,level, filters, order, alpha); 
        [en, S] = Wammerstein_adapt(un, dn, S);  
        coeffs.wammer = S.coeffs;

        err_sqr = en.^2;

        fprintf('Run time = %.2f s \n',toc);

        % Plot MSE       
        figure(fig);
        subplot(211);
        q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
        hold on; grid on;
        plot((1:iter)/1024,10*log10(MSE), 'DisplayName', ['WAMM Or.', num2str(order)]);
        axis([0 iter/1024 -40 5]);
        legend('show');

        NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
        fprintf('NMSE = %.2f dB\n', NMSE);

        % Plot NMSE
        subplot(212);
        hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['WAMM Or.', num2str(order), ': ', num2str(num2str(NMSE, '%0.2f'))] );
        axis([0 iter/1024 -20 10]);
        legend('show');

        fprintf('\n');

        %% HAMMERSTEIN
        fprintf('HAMMERSTEIN\n');
        fprintf('-------------------------------------------------------------\n');
        fprintf('Order = %d, sys_len = %d \n', order, M_hamm);        
        fprintf('mu_p = %.1f, mu_w = %.1f, alpha = %.2f \n', mu_p, mu_w,alpha)

        tic;
        S = Hammerstein_NLMS_init(order, M_hamm, [mu_p mu_w] ,leak, alpha); 
        [en, S] = Hammerstein_NLMS_adapt(un, dn, S);  
        coeffs.hammer = S.coeffs;

        err_sqr = en.^2;

        fprintf('Run time = %.2f s \n',toc);

        % Plot MSE       
        figure(fig);
        subplot(211);
        q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
        hold on; grid on;
        plot((1:iter)/1024,10*log10(MSE), 'DisplayName', ['HAMM Or.', num2str(order)]);
        axis([0 iter/1024 -40 5]);
        legend('show');

        NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
        fprintf('NMSE = %.2f dB\n', NMSE);

        % Plot NMSE
        subplot(212);
        hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['HAMM Or.' num2str(order),': ', num2str(num2str(NMSE, '%0.2f'))] );
        axis([0 iter/1024 -20 10]);
        legend('show');

        fprintf('\n');
        
        
        %% Adding title and labels to plots
        gainz = input_gain;
        figure(fig);
        subplot(211);
        title(['MSE (input type: ', input_type{i}, ', gain: ', gainz{j}, '%)']);
        ylabel('MSE dB'); grid on;
        legend('show');
        axis([0 iter/1024 -100 25]);

        subplot(212);
        title(['NMSE dB']);
        axis([0 iter/1024 -30 10]);
        xlabel('Number of iterations (\times 1024 input samples)'); 
        ylabel('Cumulative NMSE'); grid on;
        legend('show');
        
%          savefig(fig,['figures/new_',input_type{i},'_',input_gain{j},'.fig'])
        
    end
end

% %% fix plots
% if (infile == 'm') | (infile == 'f')
%     figure(MSEfig);
%     axis([0 iter/1024 -80 25]);
%     
%     figure(NMSEfig);
%     axis([0 iter/1024 -40 20]);
% end