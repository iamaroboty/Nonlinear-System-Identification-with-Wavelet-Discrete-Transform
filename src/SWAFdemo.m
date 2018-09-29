% SWAFdemo          Subband Wavelet-domain Adaptive Filter Demo
% 
% by A. Castellani & S. Cornell [Universit� Politecnica delle Marche]

addpath 'Common';             % Functions in Common folder
clear all; 
% close all;

% Adaptive filter parameters

mu = 0.1;                      % Step size
M = 256;                         % Length of unknown system response
<<<<<<< HEAD
level = 8;                       % Levels of Wavelet decomposition
=======
level = 1;                       % Levels of Wavelet decomposition
>>>>>>> 336af558b6adeae35bcd029df61ef6e6d9ac60db
wtype = 'db1';                   % Wavelet family

% Run parameters
iter = 1.0*80000;                % Number of iterations
b = load('h1.dat');              % Unknown system (select h1 or h2)
b = b(1:M);                      % Truncate to length M

%TEST: unknown system as just delay of "a" samples. Only with a integer multipler of 4, this works properly.
<<<<<<< HEAD
a = 1;     
b =[zeros(a,1); 1; zeros(M-a-1,1)];
=======
a = 0;     
b =[zeros(a,1); 1;1;1; zeros(M-a-4,1); 1];
>>>>>>> 336af558b6adeae35bcd029df61ef6e6d9ac60db

tic;

% Adaptation process

fprintf('Wavelet type: %s, levels: %d, step size = %f \n', wtype, level, mu);
[un,dn,vn] = GenerateResponses(iter,b,sum(100*clock),1,40); %iter, b, seed, ARtype, SNR
S = SWAFinit(M, mu, level, wtype);     % Initialization
S.unknownsys = b; 
[en, S] = SWAFadapt(un, dn, S);                 % Perform WSAF Algorithm

err_sqr = en.^2;
    
fprintf('Total time = %.3f mins \n',toc/60);
<<<<<<< HEAD
% dwtmode('zpd')
% [B, L] = wavedec(b, level, wtype);
% % W = WaveletMat_nL(M, level, wtype);
% % B = W*b;
% cD = detcoef(B, S.L, 'cell');
% cA = appcoef(B, S.L, wtype);
% 
% % % Plot system ID differencies
% if level == 2
%     subplot(level +1, 1, 1);
%     stem([cD{1}, S.coeffs{1}]); legend('Actual','Estimated');
%     title('cD1'); grid on;
%     subplot(level +1, 1, 2);
%     stem([cD{2}, S.coeffs{2}(:,2)]); legend('Actual','Estimated');
%     title('cD2'); grid on;
%     subplot(level +1, 1, 3);
%     stem([cA, S.coeffs{2}(:,1)]); legend('Actual','Estimated');
%     title('cA'); grid on;
%     ax=axes('Units','Normal','Position',[.1 .1 .85 .85],'Visible','off');
%     set(get(ax,'Title'),'Visible','on')
%     title('System identification. Trasformed domain');
%     h=get(ax,'Title');
% 
%     b_est = waverec([S.coeffs{2}(:,1); S.coeffs{2}(:,2); S.coeffs{1}], S.L, wtype);
%     figure;
%     stem([b, b_est]);
%     legend('Actual','Estimated');
%     title('System identification. Time domain');grid on;
% 
% elseif level == 1
%     subplot(2, 1, 1);
%     stem([cD{1}, S.coeffs{1}(:,2)]); legend('Actual','Estimated');
%     title('cD'); grid on;
%     subplot(2, 1, 2);
%     stem([cA, S.coeffs{1}(:,1)]); legend('Actual','Estimated');
%     title('cA'); grid on;
%     ax=axes('Units','Normal','Position',[.1 .1 .85 .85],'Visible','off');
%     set(get(ax,'Title'),'Visible','on')
%     title('System identification. Trasformed domain');
%     h=get(ax,'Title');
% 
%     b_est = waverec([S.coeffs{1}(:,1); S.coeffs{1}(:,2)], S.L, wtype);
%     figure;
%     stem([b, b_est]);
%     legend('Actual','Estimated');
%     title('System identification. Time domain');grid on;
% else
%     fprint('Set Level either 1 or 2\n');
% end
=======
dwtmode('sym')
[B, L] = wavedec(b, level, wtype);
% W = WaveletMat_nL(M, level, wtype);
% B = W*b;
cD = detcoef(B, S.L, 'cell');
cA = appcoef(B, S.L, wtype);


%% time domain parameters
fs = length(un)/2; % samples per sec
freq = fs/1000; % frequency
dt = 1/fs; 

%% impulse response
delta = [0; 1; zeros(fs-2,1)];
figure; 
subplot(1,2,1)
stem(delta);
title('Input Signal'); 
axis([0 10 -1.5 1.5])
out_resp = SWAtest(delta, S); 
subplot(1,2,2)
stem(out_resp);
title('Output Signal-Estimated System vs True');
hold on; 
true = conv(delta, b);
stem(true); 
axis([0 512 -1.5 1.5])



%% sine test signal 

amplitude = 1; 
padlength = 10;
input_sine = amplitude*sin(2*pi*freq*(0:dt:1-padlength*dt));
input_sine = padarray(input_sine, [0,padlength/2], 'pre'); 
input_sine = padarray(input_sine, [0,padlength/2], 'post'); 

figure; 
subplot(2,2,1)
plot(input_sine);
title('Input Signal'); 
out_sine = SWAtest(input_sine, S); 
subplot(2,2,2)
plot(out_sine);
title('Output Signal-Estimated System vs True');
hold on; 
true = conv(input_sine, b);
plot(true); 

%% FFT 
N = 2*fs;

faxis = linspace(-fs/2,fs/2,N);


subplot(2, 2, 3);
fft_true = abs(fft(input_sine, N)/N);
plot(faxis, fftshift(fft_true)); 
xlabel('Frequency');

subplot(2, 2, 4);
fft_out_est = abs(fft(out_sine, N)/N);
plot(faxis, fftshift(fft_out_est)); 
xlabel('Frequency');
hold on; 
fft_out_true = abs(fft(true,N)/N);
plot(faxis, fftshift(fft_out_true));





figure;
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
>>>>>>> 336af558b6adeae35bcd029df61ef6e6d9ac60db


figure;                          % Plot MSE
q = 0.99; MSE = filter((1-q),[1 -q],err_sqr);
hold on; plot((0:length(MSE)-1)/1024,10*log10(MSE));
axis([0 iter/1024 -60 10]);
xlabel('Number of iterations (\times 1024 input samples)'); 
ylabel('Mean-square error (with delay)'); grid on;
