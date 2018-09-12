function [en,S] = WSAFadapt(un,dn,S)

% WSAFadapt         Wavelet Subband Adaptive Filter (WSAF)
%
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal

M = length(S.coeffs);
mu = S.step;                        % Step Size
alpha = S.alpha;                 % Small constant
AdaptStart = S.AdaptStart;
W = S.W;                        % Transform Matrix
L = S.levels;
w = S.coeffs;                       % Adaptive filtering

x = zeros(M, 

ITER = length(un);
yn = zeros(1,ITER);                 % Initialize output sequence to zero
en = zeros(1,ITER);                 % Initialize error sequence to zero

if isfield(S,'unknownsys')
    b = S.unknownsys;
    norm_Wb = norm(W*b);
    eml = zeros(1,ITER);
    ComputeEML = 1;    
else
    ComputeEML = 0;
end

for n=1:ITER
    x = [un(n); x(1:end-1)];  % Fullband input vector for band partitioning
    y = [dn(n); y(1:end-1)];  % Fullband desired response vector for band partitioning
    a = [f(n); a(1:end-1)];  
    if mod(n,2)== 0         %1 level         
        U = [a'*H; U(1:end-1,:)]; %U(:,1) = cA, U(:,2) = cD
    if mod(n,4) == 0
        V = [U(1:2)*H; V(1:end-1,:)]; %V(:,1) = cA2, V(:,2) = cD2
    end
    end
end



S.coeffs = w;
if ComputeEML == 1
    S.eml = eml;
%     S.UDerr = UDerr;
end


    
