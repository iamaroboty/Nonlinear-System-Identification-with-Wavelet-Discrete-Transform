%% DISCRETE WAVELET DECOMPOSITION TESTING POLYPHASE

% Testing Signal
% clear all
close all

d = 256;        %Total signal length
t=0:0.001:10;
xn=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
xn = xn(1:d)';
% f=f(1:256);
%f = [1; -10; 324; 48; -483; 4; 7; -5532; 34; 8889; -57; 54];
%d=length(f);

% d = 512;
% f = load('h1.dat');         % Unknown system (select h1 or h2)
% f = f(1:d);                 % Truncate to length M

% f = zeros(d,1);
% f(3) = 1;

% f = chirp((0:0.001:2),0,1,250);
% f = f(1:d);
wtype = 'db2';
level = 2;
wdt_mod = 'zpd';
U = 1;

%% Generazione dei coefficienti del filtro
[low_d,high_d,low_r,high_r] = wfilters(wtype);
W = WaveletMat_nL(d, level, low_d, high_d); % DWT transform matrix
H = [low_d', high_d'];  % filter matrix analysis
F = [low_r', high_r'];  % filter matrix synthesis
[len, ~] = size(H);     % wavelet filter length

%% QMF filter
% db = dbwavf(wtype);
% qmfdb = qmf(db); 
% H = [qmfdb; db]'*sqrt(2); 
% F = flip(H);

%% Poly DWT
N = 2;
ITER = length(xn);
[len, ~] = size(H);

x = zeros(N,1);                       % Input signal buffer
X = zeros(N,len/N);                     % Polyphase buffers (for analysis filters)
Y = zeros(len/N,N);                     % Polyphase buffers (for synthesis filters)
yn = zeros(fix(ITER/N)*N,1);

for i = 1:N
    E(:,:,i) = reshape(H(:,i),N,len/N);
        % Type-I polyphase decomposition. Polyphase components of the
        %   ith filter is represented as rows in the matrix
    R(:,:,i) = flipud(reshape(F(:,i),N,len/N));
        % Type-II polyphase decomposition. Type-I Polyphase decomposition
        %   followed by a permutation 
end

for n = 1:ITER
    x = [xn(n); x(1:end-1)];
    if mod(n,N) == 0
        
           % Analysis section

        X = [x, X(:,1:end-1)];
        for i = 1:N
            xD_ply(1,i) = sum(sum(E(:,:,i).*X));
        end

           % Synthesis section

        Y = [xD_ply; Y(1:end-1,:)];
        y = 0;
        for i = 1:N
            y = y + R(:,:,i)*Y(:,i);
        end
        yn(n:-1:n-N+1) = y;
    end
end

figure;
subplot(2,1,1); plot(0:length(xn)-1,xn);
subplot(2,1,2); plot(0:length(yn)-1,yn);