% SWAFdemo          Subband Wavelet-domain Adaptive Filter Demo
% 
% by A. Castellani & S. Cornell [Università Politecnica delle Marche]

addpath 'Common';             % Functions in Common folder
clear all; close all;

% Adaptive filter parameters

mu = 0.001;                      % Step size
M = 256;                         % Length of unknown system response
level = 1;                       % Levels of Wavelet decomposition
wtype = 'db1';                   % Wavelet family

% Run parameters
iter = 1.0*80000;                % Number of iterations
b = load('h1.dat');              % Unknown system (select h1 or h2)
b = b(1:M);                      % Truncate to length M

%TEST: unknown system as just delay of "a" samples. Only with a integer multipler of 4, this works properly.
a = 2;     
b =[zeros(a,1); 1; zeros(M-a-1,1)];

tic;

% Adaptation process

fprintf('Wavelet type: %s, levels: %d, step size = %f \n', wtype, level, mu);
[un,dn,vn] = GenerateResponses(iter,b,sum(100*clock),1,60); %iter, b, seed, ARtype, SNR
S = SWAFinit(M, mu, level, wtype);     % Initialization
S.unknownsys = b; 
[en, S] = SWAFadapt(un, dn, S);                 % Perform WSAF Algorithm

err_sqr = en.^2;
    
fprintf('Total time = %.3f mins \n',toc/60);
% dwtmode('per')
% [B, L] = wavedec(b, level, wtype);
W = WaveletMat_nL(M, level, wtype);
B = W*b;
cD = detcoef(B, S.L, 'cell');
cA = appcoef(B, S.L, wtype);

% % Plot system ID differencies
if level == 2
    subplot(level +1, 1, 1);
    stem([cD{1}, S.coeffs{1}]); legend('Actual','Estimated');
    title('cD1'); grid on;
    subplot(level +1, 1, 2);
    stem([cD{2}, S.coeffs{2}(:,2)]); legend('Actual','Estimated');
    title('cD2'); grid on;
    subplot(level +1, 1, 3);
    stem([cA, S.coeffs{2}(:,1)]); legend('Actual','Estimated');
    title('cA'); grid on;
    ax=axes('Units','Normal','Position',[.1 .1 .85 .85],'Visible','off');
    set(get(ax,'Title'),'Visible','on')
    title('System identification. Trasformed domain');
    h=get(ax,'Title');

    b_est = waverec([S.coeffs{2}(:,1); S.coeffs{2}(:,2); S.coeffs{1}], S.L, wtype);
    figure;
    stem([b, b_est]);
    legend('Actual','Estimated');
    title('System identification. Time domain');grid on;

elseif level == 1
    subplot(2, 1, 1);
    stem([cD{1}, S.coeffs{1}(:,2)]); legend('Actual','Estimated');
    title('cD'); grid on;
    subplot(2, 1, 2);
    stem([cA, S.coeffs{1}(:,1)]); legend('Actual','Estimated');
    title('cA'); grid on;
    ax=axes('Units','Normal','Position',[.1 .1 .85 .85],'Visible','off');
    set(get(ax,'Title'),'Visible','on')
    title('System identification. Trasformed domain');
    h=get(ax,'Title');

    b_est = waverec([S.coeffs{1}(:,1); S.coeffs{1}(:,2)], S.L, wtype);
    figure;
    stem([b, b_est]);
    legend('Actual','Estimated');
    title('System identification. Time domain');grid on;
else
    fprint('Set Level either 1 or 2\n');
end


figure;                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
