function [en,S] = MWSAFadapt_v2(un,dn,S)
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

for i = 1:level
    A{i} = zeros(2^(i-1)*len,2^i);
    U{i} = zeros(M,2);
    D{i} = zeros(M,2);
    eDr{i} = zeros(len,1);
    delays(i) = 2^i-1; 
    eD{i} = zeros(L(end-i),1);
end
eD{i} = zeros(1,2^level);
a = zeros(len,1); 
d = zeros(len,1);   
z = zeros(len,1);
w = zeros(M,1);
% w(1) = 1;

 t=0:0.001:1;
 un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  
 dn = un;
 tot_delay = (2^level - 1)*(len-1) +1 ;


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
    A{1} = [a, A{1}(:,1:end-1)];                   % Update buffer
    
    if ComputeEML == 1
        eml(n) = norm(b-w)/norm_b;                 % System error norm (normalized)
        u = [un(n); u(1:end-1)];
        UDerr(n) = (b-w)'*u;                       % Undisturbed error 
    end
    
    % Analysis
    Dtmp = d; 
    for i = 1:level
        if (mod(n,2^i) == 0)
            Hi = upsample(H,2^(i-1));
            U1 = (Hi'*A{i})';                              % Partitioning u(n) 
            U2 = U{i}(1:end-2^i,:);
            U{i} = [U1', U2']';                           % Subband data matrix (cA || cD)            
            dD = (H'*Dtmp)';
            D{i} = [dD; D{i}(1:end-1,:)];   %cA ; cD              
            
            if i == level
                eD{i} = dD' - U{i}'*w;
                
                if n >= AdaptStart
                    w = w + U{i}*(eD{i}./(sum(U{i}.*U{i})+alpha)')*mu; % Tap-weight adaptation
%                     S.iter = S.iter + 1;
                end
            
            else
                Dtmp = D{i}(1:len,1);               % Update buffer
                tmpcat = [];                    % Improve speed here!!
                for j = 1:2^(i+1)
                    tmpcat = cat(2,tmpcat,U{i}(j:len*2^i+(j-1),1));
                end
                A{i+1} = tmpcat;
%                 A{i+1} = [U{i}(1:len*(i+1),1), U{i}(2:len*(i+1)+1,1), U{i}(3:len*(i+1)+2,1), U{i}(4:len*(i+1)+3,1)];  %Update buffer for next level, filled with cA
                eD{i} = [eD{i}(2:end); dD(2) - U{i}(:,2)'*w];
                
                if n >= AdaptStart
                    w = w + U{i}(:,2)*(eD{i}(end)./(sum(U{i}(:,2).*U{i}(:,2))+alpha)')*mu; % Tap-weight adaptation
%                     S.iter = S.iter + 1;
                end
                
            end
        end
    end

    % Synth
    for i = level:-1:1
        if (mod(n, 2^i) == 0)
            if i == level
                eDr{i} = F*eD{i} + eDr{i};
            else
                eDr{i} = F*[eDr{i+1}(1); eD{i}(end-(len-1)*delays(end-i))] + eDr{i};
                eDr{i+1} = [eDr{i+1}(2:end); 0];
            end
        end
    end
    en(n) = eDr{i}(1);
    eDr{i} = [eDr{i}(2:end); 0];  
end

S.coeffs = w;
if ComputeEML == 1
    S.eml = eml;
    S.UDerr = UDerr;
end


