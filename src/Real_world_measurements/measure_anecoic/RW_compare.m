% Real World system identification from One Plus Two device
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]


addpath '../Common';             % Functions in Common folder

clear all;  
close all;
% rng('default');

%% Structure Hyperparameters
% Filters length
M1 = 256;
M2 = 32;

% Universal stepsize
mu = [0.1 0.1];
mu_fb = [0.1 0.5];

% For Wavelet Transform
level = 3;
filters = 'db16';

% Hammerstein and Wammerstein
% p , w
leak = [0 0];
order = 2;
M_hamm = M1 + M2;
mu_p = 0.3;
mu_w = 0.07;
alpha = 10.^-1;

% For MSAF and SAF
N = 2^level;                % Number of subbands                     
D = N/2;                    % Decimation factor for 2x oversampling
L = 8*N;                    % Length of analysis filters, M=2KN

% Iter count
iter = 0.5*80000;

input_type = {'simulated'};
input_gain = {'0'};

n = 0;
for i= 1:numel(input_type)
    for j= 1:numel(input_gain)
        n = n+1;
        fprintf('STARTING RUN (%d)/(%d), file: %s \n', n , numel(input_type)*numel(input_gain), [input_type{i}, '_', input_gain{j}]);
        %% Input Data
%         % Anecoic measurment
%         load(['univpm_ref_', input_type{i}, '.mat'])
%         un = un(:)';    % Make it row vector
%         load(['univpm_', input_type{i}, '_', input_gain{j},'.mat'])
%         dn = dn(:)';    % Make it row vector
% 
%         delay = 1024 + 100;
%         dn = dn(delay:end); 
%         un = un(1:end-delay+1); 
% 
%         dn = dn(1:iter);
%         un = un(1:iter);
% 
%         fs = 44100;
% 
%         % % Resample the data to new fsn
% %         fs_new = 8000;
% %         [P, Q] = rat(fs_new/fs);
% %         dn = resample(dn, P, Q);
% %         un = resample(un, P, Q);
% %         iter = length(un);
%         
% 
% %         Normalization of u(n) and d(n)
% %         un = un/std(dn);
% %         dn = dn/std(dn);
% 
% %         Normalization of speech signal
%         a = abs(max(dn))/sqrt(2);
%         un = un/a; dn = dn/a;

        %% Create and plot kernel simulated
        % Create Kernel mode
        kermode = 'simulated';              %modes: "delta", "randomdiag", "random", "lowpass", "simulated"

        % Create Kernel parameters
        deltapos = [5, 3];
        Ndiag = 1;
        normfreq = [0.6, 0.2];
        h1h2 = ["h1.dat", "h2.dat"];
        param = {deltapos, Ndiag, normfreq, h1h2};

        [ker1, ker2] = create_kernel(M1, M2, kermode, param);

        NL_system.M = [M1, M2];
        gains = [1 0.1]; 
        NL_system.Responses = {gains(1).*ker1, gains(2).*ker2};

        % Plotting kernel
        figure('Name', 'Simulated Kernel');
        kernel_plot(NL_system.Responses);
        [un,dn,vn] = GenerateResponses_Volterra(iter, NL_system ,sum(100*clock),4,40); %iter, b, seed, ARtype, SNR
        iter = length(un);
        % b = load('h1.dat');
        % b = b(1:M);
        % [un,dn,vn] = GenerateResponses(iter,b,sum(100*clock),1,40); %iter, b, seed, ARtype, SNR
        % [un,dn,vn] = GenerateResponses_speech(b,'SpeechSample.mat');
        
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
        fprintf('step size = %s \n', sprintf('%.2f ', mu));
        
        tic;
        S = Volterra_Init([M1, M2], mu, level, filters); 
        [en, S] = Volterra_2ord_adapt_v3(un, dn, S);
        
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
        % fprintf('MSAFTERRA\n');
        % fprintf('--------------------------------------------------------------------\n');
        % fprintf('Number of subbands, N = %d, step size = %s \n', N, sprintf('%.2f ', mu));
        % 
        % S = MSAFTERRA_Init([M1 M2],mu,N,L);
        % tic;
        % [en,S] = MSAFTERRA_adapt(un,dn,S);
        % err_sqr = en.^2;
        % 
        % fprintf('Total time = %.2f s \n',toc);
        % 
        % % Plot MSE       
        % figure(MSEfig);
        % q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
        % hold on; grid on;
        % plot((1:iter)/1024,10*log10(MSE), 'DisplayName', ['MSAFTERRA, N: ', num2str(N), ' subband']);
        % axis([0 iter/1024 -40 5]);
        % ylabel('Mean-square error'); grid on;
        % xlabel('Number of iterations (\times 1024 input samples)');
        % legend('show');
        % 
        % NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
        % fprintf('NMSE = %.2f dB\n', NMSE);
        % 
        % % Plot NMSE
        % figure(NMSEfig);
        % hold on; plot((0:iter-1)/1024, 10*log10(cumsum(err_sqr)./(cumsum(dn.^2))), 'DisplayName', ['MSAFTERRA, N: ', num2str(N), ' subband'] );
        % axis([0 iter/1024 -20 10]);
        % xlabel('Number of iterations (\times 1024 input samples)'); 
        % ylabel('Cumulative Normalized Mean-square error'); grid on;
        % legend('show');
        % 
        % fprintf('\n');


        %% SAFTERRA
        % fprintf('SAFTERRA\n');
        % fprintf('--------------------------------------------------------------------\n');
        % fprintf('Number of subbands, N = %d, Decimation factor, D = %d, step size = %s \n', N, D, sprintf('%.2f ', mu));
        % 
        % S = SAFTERRA_Init([M1 M2],mu,N,D,L);
        % tic;
        % [en,S] = SAFTERRA_adapt(un,dn,S);
        % err_sqr = en.^2;
        % 
        % fprintf('Total time = %.2f s \n',toc);
        % 
        % figure(MSEfig);
        % q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
        % hold on; plot((1:iter)/1024,10*log10(MSE), 'DisplayName', ['SAFTERRA, N: ', num2str(N), ' subband'] );
        % ylabel('Mean-square error'); grid on;
        % legend('show');
        % axis([0 iter/1024 -40 5]);
        % xlabel('Number of iterations (\times 1024 input samples)');
        % 
        % NMSE = 10*log10(sum(err_sqr)/sum(dn.^2));
        % fprintf('NMSE = %.2f dB\n', NMSE);
        % 
        % % Plot NMSE
        % figure(NMSE_fig);  
        % [NMSE, indx] = NMSE_compute(dn, en, nmse_n_points);
        % hold on; grid on;
        % plot(indx, NMSE,'DisplayName', ['SAFTERRA, N: ', num2str(N), ' subband']);
        % axis([indx(1) indx(end) -10 0]);
        % xlabel('Iteration number'); 
        % ylabel('Cumulative Normalized Mean-square error'); 
        % legend('show');
        % 
        % fprintf('\n');

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
        gainz = {'0 dB'};
        figure(fig);
        subplot(211);
        title(['MSE (input type: ', input_type{i}, ', gain: ', gainz{j}, ')']);
        ylabel('MSE dB'); grid on;
        legend('show');
        axis([0 iter/1024 -100 20]);

        subplot(212);
        title(['NMSE dB']);
        axis([0 iter/1024 -40 15]);
        xlabel('Number of iterations (\times 1024 input samples)'); 
        ylabel('Cumulative NMSE'); grid on;
        legend('show');
        
%          savefig(fig,['figures/nonlin_',input_type{i},'_',input_gain{j},'.fig'])
        
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