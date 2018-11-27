clear, clc, close all

% signal parameters
fs = 44100;         % sampling frequency, Hz
T = 5;              % signal duration, s
N = round(fs*T);    % number of samples

% noise generation
xred = rednoise(1, N);         % red (Brownian) noise
xpink = pinknoise(1, N);       % pink noise
xblue = bluenoise(1, N);       % blue noise
xviolet = violetnoise(1, N);   % violet noise

% sound presentation
soundsc(xred, fs)
pause(T+1)

soundsc(xpink, fs)
pause(T+1)

soundsc(xblue, fs)
pause(T+1)

soundsc(xviolet, fs)