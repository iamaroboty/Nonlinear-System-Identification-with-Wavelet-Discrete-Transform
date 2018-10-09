close all; clear all;
fs = 100;               % sampling frequency
t = 0:(1/fs):(10-1/fs); % time vector
S = cos(2*pi*10*t);
n = length(S);
X = fft(S);
f = (0:n-1)*(fs/n);     %frequency range
power = abs(X).^2/n;    %power
figure; subplot(2,1,1); plot(f,power)

Y = fftshift(X);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power
subplot(2,1,2); plot(fshift,powershift)

%% Wavelet QMF
% Load scaling filter associated with an orthogonal wavelet. 
wtype = 'db10';
[low_d,high_d,low_r,high_r] = wfilters(wtype);

db = dbwavf(wtype);
subplot(321); stem(db); title(sprintf('%s low-pass filter', wtype));

% Compute the quadrature mirror filter. 
qmfdb = qmf(db); 
subplot(322); stem(qmfdb); title(sprintf('QMF %s filter', wtype));

% Check for frequency condition (necessary for orthogonality):
% abs(fft(filter))^2 + abs(fft(qmf(filter))^2 = 1 at each 
% frequency. 
m = fft(db); 
mt = fft(qmfdb); 
freq = [1:length(db)]/length(db); 
subplot(323); plot(freq,abs(m)); 
title(sprintf('Transfer modulus of %s', wtype))
subplot(324); plot(freq,abs(mt)); 
title(sprintf('Transfer modulus of QMF %s', wtype))
subplot(325); plot(freq,abs(m).^2 + abs(mt).^2); 
title(sprintf('Check QMF condition for %s and QMF %s', wtype, wtype)) 
xlabel(sprintf(' abs(fft(%s))^2 + abs(fft(qmf(%s))^2 = 1', wtype, wtype))
axis([0 1 0 2])

% Editing some graphical properties,
% the following figure is generated.

% Check for orthonormality. 
df = [db;qmfdb]*sqrt(2); 
id = df*df'