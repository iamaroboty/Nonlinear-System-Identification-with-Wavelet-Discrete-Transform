clear all; close all
% Set initial signal and get filters.
x = sin(0.3*[1:451]); w = 'db9'; 
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(w);
% In fact using a slightly redundant scheme, any signal
% extension strategy works well. 
% For example use random padding.

lx = length(x); lf = length(Lo_D);
ex = [randn(1,lf) x randn(1,lf)];
% ex = [zeros(1,lf) x zeros(1,lf)];
axis([1 lx+2*lf -2 3])
subplot(211), plot(lf+1:lf+lx,x), title('Original signal')
axis([1 lx+2*lf -2 3])
subplot(212), plot(ex), title('Extended signal')
axis([1 lx+2*lf -2 3])

% Decomposition.
la = floor((lx+lf-1)/2);
ar = wkeep(dyaddown(conv(ex,Lo_D)),la);
dr = wkeep(dyaddown(conv(ex,Hi_D)),la);
% Reconstruction.
xr = idwt(ar,dr,w,lx);

% Check perfect reconstruction.
err0 = max(abs(x-xr))

%% 
level = 2;
dwtmode('zpd')
[C, L] = wavedec(x, level, w);
xr1 = waverec(C, L, w);
err1 = max(abs(x-xr1))