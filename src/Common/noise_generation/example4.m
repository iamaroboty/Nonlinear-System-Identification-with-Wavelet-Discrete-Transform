clear, clc, close all

% sampling frequency, Hz
fs = 44100;         

% noise generation
xred = rednoise(1000, 1000);         % red (Brownian) noise 
xpink = pinknoise(1000, 1000);       % pink noise
xblue = bluenoise(1000, 1000);       % blue noise
xviolet = violetnoise(1000, 1000);   % violet noise

% visualize the signals
figure(1)
colormap gray

subplot(2, 2, 1)
imagesc(xred)
title('Red noise')

subplot(2, 2, 2)
imagesc(xpink)
title('Pink noise')

subplot(2, 2, 3)
imagesc(xblue)
title('Blue noise')

subplot(2, 2, 4)
imagesc(xviolet)
title('Violet noise')