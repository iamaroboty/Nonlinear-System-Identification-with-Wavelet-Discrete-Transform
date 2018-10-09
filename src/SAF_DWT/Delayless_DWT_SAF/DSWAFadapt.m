function [en,S] = DSWAFadapt(un,dn,S)
% SWAFadapt         Wavelet-transformed Subband Adaptive Filter (WAF)                 
%
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal

M = S.length;                     % Unknown system length (Fullband adaptive filter)
mu = S.step;                      % Step Size
AdaptStart = S.AdaptStart;        % Transient
alpha = S.alpha;                  % Small constant (1e-6)
H = S.analysis;                   % Analysis filter bank
F = S.synthesis;                  % Synthesis filter bank
[len, ~] = size(H);               % Wavelet filter length
level = S.levels;                 % Wavelet Levels
L = S.L;                          % Wavelet decomposition Length, subfilter length [cAn cDn cDn-1 ... cD1 M]
UpdateRate = S.UpdateRate;
Hadam = hadamard(2^level);

dif = L(1) - M/2;

% E = S.Polyphase;


% Init Arrays
for i= 1:level
    U.cD{i} = zeros(L(end-i),1);    
    U.cA{i} = zeros(L(end-i),1);    
    w{i} = zeros(L(end-i),1);       % Subband adaptive filter coefficient, initialize to zeros    
end 
w{i} = zeros(L(end-i),2);           % Last level has 2 columns, cD and cA
eD{i} = zeros(2,1);                 % Last level has 2 columns, cD and cA
W = zeros(1,M);                     % Fullband adaptive filter
U.tmp = zeros(len,1);
U.Z = zeros(2,1);

x = zeros(M,1);                    % Fullband tap-input vector
u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  
e = zeros(len,1);                 % Tapped-delay line of error

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


for n = 1:ITER    
    x = [un(n); x(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    
    u = [un(n); u(1:end-1)];
    
    en(n) = dn(n) - W*x;         % Error signal
    e = [en(n); e(1:end-1)];
          
    % Analysis Bank
    U.tmp = u;
    for i = 1:level
        if mod(n,2^i) == 0
            U.Z = H'*U.tmp;
            U.cD{i} = [U.Z(2); U.cD{i}(1:end-1)]; 
            U.cA{i} = [U.Z(1); U.cA{i}(1:end-1)];
            U.tmp = U.cA{i}(1:len);                        
            
            if i == level
                eD{i} = H'*e;
                
                if n >= AdaptStart(i)
                    w{i} = w{i} + mu*[U.cA{i},U.cD{i}].*eD{i}'./(sum([U.cA{i},U.cD{i}].*[U.cA{i},U.cD{i}])+alpha); 
                end 
            else
%                 eD{i} = [eD{i}(2:end); Y.cD{i}(1) - U.cD{i}'*w{i}]; 
%                 
%                 if n >= AdaptStart(i)
% %                     pwr{i} = beta(i)*pwr{i}+ (1-beta(i))*(U.cD{i}.*U.cD{i});
% %                     w{i} = w{i} + (mu*eD{i}(end)/((sum(pwr{i}) + alpha)))*U.cD{i};
%                     w{i} = w{i} + (mu*eD{i}(end)/(U.cD{i}'*U.cD{i} + alpha))*U.cD{i};
%                 end
            end           
            S.iter{i} = S.iter{i} + 1;                
        end
    end    
    
    if (n >= AdaptStart+M) && (mod(n,UpdateRate)==0)          
%         W = (waverec([w{i}(:,1); w{i}(:,2)], L, F(:,1), F(:,2)))';
        W = 1/2.* Hadam * [w{i}(1:M/2,2), w{i}(1:M/2,1)]';
        W = W(:)';
        
    end
   
end

S.FULLcoeffs = W;
S.SUBcoeffs = w;
end

