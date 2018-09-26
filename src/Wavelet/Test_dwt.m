%% DISCRETE WAVELET DECOMPOSITION TESTING

% Testing Signal
clear all; close all

d = 256;        %Total signal length
t=0:0.001:10;
f=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
f = f(1:d)';
% f=f(1:256);
%f = [1; -10; 324; 48; -483; 4; 7; -5532; 34; 8889; -57; 54];
%d=length(f);

% d = 512;
% f = load('h1.dat');         % Unknown system (select h1 or h2)
% f = f(1:d);                 % Truncate to length M

% f = zeros(d, 1);
% f(129) = 1;

wtype = 'db1';
level = 1;

%% Generazione dei coefficienti del filtro
[low_d,high_d,low_r,high_r] = wfilters(wtype);
W = WaveletMat_nL(d, level, low_d, high_d); % DWT transform matrix
H = [low_d', high_d'];  % filter matrix analysis
F = [low_r', high_r'];  % filter matrix synthesis
[len, ~] = size(H);     % wavelet filter length

%% Decomposizone con funzione matlab
dwtmode('zpd')
tic
[C, L] = wavedec(f, level, wtype);
fr = waverec(C, L, wtype);
toc
D = detcoef(C,L,'cells');
A = appcoef(C, L, wtype);
errC = max(abs(f-fr))

%% Decomposizione con Matrice W
tic
% Zeropads the matrix border
padsize = len/2-1;
W_pad = padarray(W, [padsize padsize], 'both');
% W_pad = [zeros(d,1), W, zeros(256,1)];
% W_pad = [zeros(d+2,1)'; W_pad];
% W_pad = [W_pad(1:d/2,:); zeros(d+2,1)'; W_pad(d/2+1:end,:)]; %this works as matlab function
f_pad = padarray(f, padsize, 'both');

Z = W_pad*f_pad;
Zr = W_pad'*Z;
toc

Zr = Zr(1+padsize:end-padsize);
errZ = max(abs(f-Zr))

for dval=1:level 
    zD{dval} = Z((d/(2^(dval)))+1:(d/(2^(dval-1)))); 
    if (dval == level) % only last level for cA
        zA{dval} = Z(1:(d/(2^(dval)))); 
    end 
end

if length(Z) == length(C)
    diff_ZC = max(abs(C-Z)) %diff from C and Z
else    
    fprintf(['-------------------------------------------\n',...
        'WARNING: C and Z have different length \n',...
        '-------------------------------------------\n']);
end

% stem([Z, C]); legend('Z','MatLab');
% title('Z-C'); grid on;

%% Decomposizione manuale
% % Analysis
% [McA,McD] = dwt(f,low_d,high_d);       %level 1
% [McA2,McD2] = dwt(McA,low_d,high_d);   %level 2
% 
% % Synthesis
% McR2 = idwt(McA2,McD2,low_r,high_r);
% McR = idwt(McR2,McD,low_r,high_r);
% err0 = max(abs(f-McR))

%% Border Effect 
lx = length(f); lf = length(low_d);
% ex = [randn(1,lf) f' randn(1,lf)];
ex = [zeros(1,lf) f' zeros(1,lf)];
% dwtmode('zpd')

% Decomposition.
la = floor((lx+lf-1)/2);
ar = wkeep(dyaddown(conv(ex,low_d)),la);
dr = wkeep(dyaddown(conv(ex,high_d)),la);
% Reconstruction.
xr = idwt(ar,dr,wtype,lx); %central portion

% Check perfect reconstruction.
err1 = max(abs(f'-xr))

%% Real time DWT
tic
delay = length(H)-1;            %total delay (analysis + synthesis)
x = zeros(length(H)^level,1);
b = zeros(length(F)^level,1);
yn = zeros(d+delay, 1);
fpad = [f; zeros(delay, 1)];

cD = zeros(d/2, 1);
cA = zeros(d/2, 1);

for n = 1:d+delay
    x = [fpad(n); x(1:end-1)];
    if mod(n,2) == 0
        xD = H'*x;
        cD = [cD(2:end); xD(2)];
        cA = [cA(2:end); xD(1)];
        b = F*xD + b;
    end
    yn(n) = b(1);
    b = [b(2:end); 0];
end
toc
err2 = max(abs(f-yn(1+delay:end)))

%% Plotting Stuff
% figure; 
% subplot(3,1,1);
% plot(f); title('Original'); axis([1 d -inf inf]);
% subplot(3,1,2);
% plot([zD{1}, D{1}, cD]); axis([1 d/2 -inf inf]);
% title('cD'); legend('Matrix W', 'MATLAB', 'Real-Time'); grid on;
% subplot(3,1,3);
% plot([zA{1}, A, cA]); axis([1 d/2 -inf inf]);
% title('cA'); legend('Matrix W', 'MATLAB', 'Real-Time'); grid on;

