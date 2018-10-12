function [en,S] = SWAFtest_WAVPACK_v2(un, S, Ovr)

% SWAFadapt         Wavelet-transformed Subband Adaptive Filter (WAF)                 
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
L = S.L.*Ovr./2;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

% Init Arrays
for i= 1:level
     
       U.c{i} = zeros(L(end-i),4^(i));        
       Y.c{i} = zeros(L(end-i),4^(i));  
       eD{i} = zeros(L(end-i),4^(i-1));      % Error signa, transformed domain
       eDr{i} = zeros(len,4^(i-1));          % Error signal, time domain
       delays(i) = 4^i-1;                    % Level delay for synthesis
    
        
end 
w = S.coeffs;           % Last level has 2 columns, cD and cA
eD{i} = zeros(1,4^i);   % Last level has 2 columns, cD and cA

pwr = w;
beta = 1./L(2:end-1);

u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  
y = zeros(len,1);                 % Tapped-delay line of desired response (Analysis FB)

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


% % ONLY FOR TESTING PURPOSE
%  t=0:0.001:1;
%  un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  
%  ITER = length(un);
% Testing freezed filters
% w{1} = zeros(L(end-1),2);
% 
%  w(1,1:end) = 1;

for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
      
    % Analysis Bank
    U.tmp = u;
     
    for i = 1:level
        if mod(n,4^i/(Ovr)) == 0
            if (i==1 && Ovr == 2)
                HH = H./sqrt(2);
            else
                HH = H;
            end
            
            U.Z = HH*U.tmp; % column U.cD U.cA
                    
            [rows, cols] = size(U.Z);
            
            indx = 1;
            
            for col=1:cols
                for row=1:rows 
                    
                U.c{i}(:,indx) = cat(1,U.Z(row,col), U.c{i}(1:end-1, indx));
                           
                indx=indx+1;
                end  
            end
            
            U.tmp = U.c{i}(1:len,:);    
          
            if i == level
                
                eD{i} = sum((U.c{i}).*w);
                       
            end           
            S.iter{i} = S.iter{i} + 1;                
        end
    end    


    % Synthesis Bank
    for i = level:-1:1
            if (i==1 && Ovr == 2)
                FF = F./sqrt(2);
            else
                FF = F;
            end
       if i == level
            if mod(n,4^i/(Ovr)) == 0
                    eDr{i} = FF*eD{i}' + eDr{i};                       
            end
        else
            if mod(n,4^i/(Ovr)) == 0     
                
                indx = 1; 
                
                for col = 1:2:2^i-1
                    eDr{i}(:,indx) = FF*eDr{i+1}(1,col:col+1)' + eDr{i}(:,indx);
                    indx = indx +1;                
                end
                eDr{i+1} = [eDr{i+1}(2:end,:); zeros(1,size(eDr{i+1},2))];
            end            
        end
    end   
    en(n) = eDr{i}(1);
    eDr{i} = [eDr{i}(2:end); 0];           
end

en = en(1:ITER);
S.coeffs = w;
end

