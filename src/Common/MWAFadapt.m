function [en,S] = MWAFadapt(un,dn,S, dwt_mode)

% WSAFadapt         Wavelet-transformed Subband Adaptive Filter (WAF)
%                   Transformation with an orthogonal matrix W.                    
%
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal

M = length(S.coeffs);
mu = S.step;                     % Step Size
beta = S.beta;                   % Forgettig factor
AdaptStart = S.AdaptStart;
W = S.W;                         % Transform Matrix
level = S.levels;                  % Wavelet Levels
wtype = S.wtype;                 % Wavelet family   
b = S.unknownsys;

L = zeros(level,1);                   % Each subband lenght [cAn cDn cDn-1 ... cD1]
H = S.analysis;
[len, ~] = size(H);

u = zeros(M,1);
y = zeros(M,1);

ITER = length(un);
en = zeros(1,ITER);                 % Initialize error sequence to zero

dwtmode(dwt_mode)

if dwt_mode == 'per'
    for i= 1:level
        L = [M/(2^i); L(1:end-1)];
    end
    L = [L(1); L; M]';
else
    L = [M; zeros(level,1)];
    for i= 1:level
        L = [floor((L(1)+len-1)/2); L(1:end-1)];
    end
    L = [L(1); L]';
end           
S.L = L;
    
% Init arrays
for i = 1:level
    w{i} = zeros(L(end-i),1); 
end
w{i} = zeros(L(end-i),2); 
eD = w;

% ONLY FOR TESTING PURPOSE
% t=0:0.001:1;
% dn=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
% [UU,LL] = wavedec(flip(dn(1:256)), level, wtype);
% UU = W*flip(un(1:256))';
for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    y = [dn(n); y(1:end-1)];        % Desired response vector        
    
    for i= 1:level
        if mod(n,2^i) == 0
            if i == 1
                % Wavelet transform
                [U, ~] = wavedec(u, level, wtype);
                [Y, ~] = wavedec(y, level, wtype);

                % Decomposition
                UcD = detcoef(U, L, 'cells');
                YcD = detcoef(Y, L, 'cells');

                UcA = appcoef(U, L, wtype);
                YcA = appcoef(Y, L, wtype);
            end
            
            if i == level
                eD{i} = [[YcA(1), YcD{i}(1)] - sum([UcA, UcD{i}].*w{i}); eD{i}(1:end-1,:)];
                
                if n >= AdaptStart
                    w{i} = w{i}+(mu*eD{i}(1,:)./(sum([UcA, UcD{i}].*[UcA, UcD{i}])+0.00001).*[UcA, UcD{i}]);
                end
                
                % Synthesis
%                ew = [ecA; ecD4; ecD3; ecD2; ecD1];  
                ew = cat(1,eD{level-1:-1:1});
                ew = [reshape(eD{level},[L(1)*2,1]); ew];
                                                                      
                z = waverec(ew, L, wtype);
                en(n-2^i+1:n) = flip(z(1:2^i));
            else
                eD{i} = [YcD{i}(1) - sum(UcD{i}.*w{i}); eD{i}(1:end-1)];
                
                if n >= AdaptStart
                    w{i} = w{i}+(mu*eD{i}(1)/(UcD{i}'*UcD{i} + 0.00001)*UcD{i});
                end
            end
        end
    end
    
%     if mod(n,5000)== 0
%         plot(10*log10(en(1:n).^2));
%         xlabel('Number of iteration'); 
%         ylabel('Live MSE error (dB)');linkdata on    %Live plotting      
%     end
    
end


en = en(1:ITER);
S.coeffs = w;
end


    
