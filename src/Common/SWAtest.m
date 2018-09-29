function [y_out] = SWAtest(un,S)
% TESTING FOR 2 LEVELS, EASY EXTENSION TO GENERAL LEVEL AND WAVELET
% FUNCTIONS
% SWAFadapt         Wavelet-transformed Subband Adaptive Filter (WAF)                 
%
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal

H = S.analysis;
F = S.synthesis;
[len, ~] = size(H);               % Wavelet filter length
level = S.levels;                 % Wavelet Levels
L = S.L;                          % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1]

for i= 1:level
    U{i} = zeros(L(end-i),2);         % Tapped-delay lines of adaptive subfilters, wavelet coefficient C{1} = cA, C{2}:C{end} = cDs
end  

w = S.coeffs;

Z = zeros(2,1);

u = zeros(len,1); 
 % Tapped-delay line of input signal (Analysis FB)  
             % Tapped-delay line of desired response (Analysis FB)

ITER = length(un);
y_out = zeros(1,ITER);

for n = 1:ITER
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    if mod(n,2) == 0
        U{1} = [u'*H; U{1}(1:end-1,:)];
        Y = sum(U{1}.*w{1});
      
        Z = F*Y';
    end
    y_out(n) = Z(1);
    Z = [Z(2:end); 0]; 
end    

y_out = y_out(1:ITER);

