function [en,S] = MWSAFadapt(un,dn,S)
% Wavelet-Decomposition Subband Adaptive Filter (WAF)                 
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal

M = S.length;                     % Unknown system length (Equivalent adpative filter lenght)
mu = S.step;                      % Step Size
AdaptStart = S.AdaptStart;        % Transient
alpha = S.alpha;                  % Small constant (1e-6)
H = S.analysis;                   % Analysis filter bank
F = S.synthesis;                  % Synthesis filter bank

[len, ~] = size(H);               % Wavelet filter length
level = S.levels;                 % Wavelet Levels
L = S.L;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

U = zeros(M,2^level);                      % Adaptive filtering
a = zeros(len,1); d = zeros(len,1); A = zeros(len,2^level);    % Analysis filtering
z = zeros(len,1);
w = zeros(M,1);

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


if isfield(S,'unknownsys')
    b = S.unknownsys;
    norm_b = norm(b);
    eml = zeros(1,ITER);
    ComputeEML = 1;
    u = zeros(M,1);
else
    ComputeEML = 0;
end
	
for n = 1:ITER
    
    d = [dn(n); d(1:end-1)];                       % Update tapped-delay line of d(n)
    a = [un(n); a(1:end-1)];                       % Update tapped-delay line of u(n)
    A = [a, A(:,1:end-1)];                         % Update buffer
    if ComputeEML == 1
        eml(n) = norm(b-w)/norm_b;                 % System error norm (normalized)
        u = [un(n); u(1:end-1)];
        UDerr(n) = (b-w)'*u;                       % Undisturbed error 
    end
    
    if (mod(n,2^level)==0)                               % Tap-weight adaptation at decimated rate
        U1 = (H'*A)';                              % Partitioning u(n) 
        U2 = U(1:end-2^level,:);
        U = [U1', U2']';                           % Subband data matrix
        dD = H'*d;                                 % Partitioning d(n) 
        eD = dD - U'*w;                            % Error estimation
        if n >= AdaptStart
            w = w + U*(eD./(sum(U.*U)+alpha)')*mu; % Tap-weight adaptation
            S.iter = S.iter + 1;
        end
        z = F*eD + z;                                       
        en(n-2^level+1:n) = z(1:2^level); 
        z = [z(2^level+1:end); zeros(2^level,1)];
    end
                          
end

S.coeffs = w;
if ComputeEML == 1
    S.eml = eml;
    S.UDerr = UDerr;
end


    


