function [en,S] = Wammerstein_adapt(un,dn,S)
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
level = S.levels;                 % Wavelet Levels

H = create_multilevel_bank(H, level);
F = flip(H);

S.analysis_bank = H;
S.synthesis_bank = F;

[len, ~] = size(H);               % Wavelet filter length
L = S.L;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

order = S.order; 
xp = zeros(order, 1);  
w = zeros(M,1);
p = zeros(order, 1); % non linearity coeffs vector 
p(1) = 1; 
X = zeros(M, order); 

U = zeros(M,2^level);                      % Adaptive filtering
a = zeros(len,1); d = zeros(len,1); A = zeros(len,2^level);    % Analysis filtering
z = zeros(len,1);

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero

	
for n = 1:ITER
    
    for i = 1:order % build the coeff vector (taylor expansion)        
         xp = [xp(2:end); un(n)^i];                   
    end     
    X = cat(1, xp', X(1:end-1,:)); 
    x_nl = X*p;    
    
    d = [dn(n); d(1:end-1)];                       % Update tapped-delay line of d(n)
    A = [x_nl(1:len), A(:,1:end-1)];                         % Update buffer
    
    if n >= AdaptStart
        p = p + (mu(1)*en(n)*X'*w)/((X'*w)'*(X'*w) + alpha);    
    end
    
    if (mod(n,2^level)==0)                               % Tap-weight adaptation at decimated rate
        U1 = (H'*A)';                              % Partitioning u(n) 
        U2 = U(1:end-2^level,:);
        U = [U1', U2']';                           % Subband data matrix
        dD = H'*d;                                 % Partitioning d(n) 
        eD = dD - U'*w;                            % Error estimation
        if n >= AdaptStart
            w = w + U*(eD./(sum(U.*U)+alpha)')*mu(2); % Tap-weight adaptation
%             S.iter = S.iter + 1;
        end
        z = F*eD + z;                                       
%         en(n-2^level+1:n) = z(1:2^level); 
%         z = [z(2^level+1:end); zeros(2^level,1)]; 
    end
    en(n) = z(1);
    z = [z(2:end); 0];                     
end

S.coeffs{1} = w;
S.coeffs{2} = p; 
end


    


