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

% H = S.analysis;
% F = S.synthesis;
% [len, ~] = size(H);               % Wavelet filter length
% level = S.levels;                 % Wavelet Levels
% L = S.L;                          % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1]
% 
% for i= 1:level
%     U{i} = zeros(L(end-i),2);         % Tapped-delay lines of adaptive subfilters, wavelet coefficient C{1} = cA, C{2}:C{end} = cDs
% end  
% 
% w = S.coeffs;
% 
% Z = zeros(2,1);
% 
% u = zeros(len,1); 
%  % Tapped-delay line of input signal (Analysis FB)  
%              % Tapped-delay line of desired response (Analysis FB)
% 
% ITER = length(un);
% y_out = zeros(1,ITER);
% 
% for n = 1:ITER
%     u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
%     if mod(n,2) == 0
%         U{1} = [u'*H; U{1}(1:end-1,:)];
%         Y = sum(U{1}.*w{1});
%       
%         Z = F*Y';
%     end
%     y_out(n) = Z(1);
%     Z = [Z(2:end); 0]; 
% end    
% 
% y_out = y_out(1:ITER);



%M = S.length;                     % Unknown system length (Equivalent adpative filter lenght)
%mu = S.step;                      % Step Size
%AdaptStart = S.AdaptStart;        % Transient
%alpha = S.alpha;                  % Small constant (1e-6)
H = S.analysis;                   % Analysis filter bank
F = S.synthesis;                  % Synthesis filter bank
[len, ~] = size(H);               % Wavelet filter length
level = S.levels;                 % Wavelet Levels
L = S.L;                          % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

% Init Arrays
for i= 1:level
    U.cD{i} = zeros(L(end-i),1);    
    U.cA{i} = zeros(L(end-i),1);    
    Y.cD{i} = zeros(L(end-i),1);
    Y.cA{i} = zeros(L(end-i),1);
    eD{i} = zeros(L(end-i),1);      % Error signa, transformed domain
    eDr{i} = zeros(len,1);          % Error signal, time domain
    delays(i) = 2^i-1;              % Level delay for synthesis
    w{i} = S.coeffs{i};%zeros(L(end-i),1);       % Subband adaptive filter coefficient, initialize to zeros
end 
w{i} = S.coeffs{i}; %zeros(L(end-i),2);           % Last level has 2 columns, cD and cA
eD{i} = zeros(1,2);                 % Last level has 2 columns, cD and cA
U.tmp = zeros(len,1);
Y.tmp = zeros(len,1);
U.Z = zeros(2,1);
Y.Z = zeros(2,1);

u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  
y = zeros(len,1);                 % Tapped-delay line of desired response (Analysis FB)

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


% % ONLY FOR TESTING PURPOSE
% t=0:0.001:1;
% un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  

% Testing freezed filters
% w{1} = zeros(L(end-1),2);
% w{1}(1,:) = 1;

for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    %y = [dn(n); y(1:end-1)];        % Desired response vector        

    % Analysis Bank
    U.tmp = u;
    for i = 1:level
        if mod(n,2^i) == 0
            U.Z = H'*U.tmp;
            U.cD{i} = [U.Z(2); U.cD{i}(1:end-1)]; 
            U.cA{i} = [U.Z(1); U.cA{i}(1:end-1)];
            U.tmp = U.cA{i}(1:len);
            
            if i == level
                eD{i} = sum(([U.cA{i}, U.cD{i}]).*w{i});

                
            else
                eD{i} = [eD{i}(2:end); U.cD{i}'*w{i}]; 

            end
           
            S.iter{i} = S.iter{i} + 1;
                
        end
    end    


    % Synthesis Bank
    for i = level:-1:1
        if i == level
            if mod(n,2^i) == 0
                eDr{i} = F*eD{i}' + eDr{i};
            end
        else
            if mod(n,2^i) == 0                
                eDr{i} = F*[eDr{i+1}(1); eD{i}(end-(len-1)*delays(end-i))] + eDr{i};
                eDr{i+1} = [eDr{i+1}(2:end); 0];
            end            
        end
    end   
    en(n) = eDr{i}(1);
    eDr{i} = [eDr{i}(2:end); 0];           
end

y_out = en(1:ITER);

end




