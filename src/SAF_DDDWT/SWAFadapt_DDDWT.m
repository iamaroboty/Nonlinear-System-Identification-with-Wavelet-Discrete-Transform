function [en,S] = SWAFadapt_DDDWT(un,dn,S)
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
level = S.levels;                 % Wavelet Levels
H = S.analysis;
F = S.synthesis;                     
L = S.L;
len = size(H,1);

% Init Arrays
for i= 1:level               
    U.cA{i} = zeros(L(end-i),1); 
    U.cD1{i} = zeros(L(end-i),1);  
    U.cD2{i} = zeros(L(end-i),1);  
    
    Y.cA{i} = zeros(L(end-i),1);
    Y.cD1{i} = zeros(L(end-i),1);
    Y.cD2{i} = zeros(L(end-i),1);    

    eD{i} = zeros(L(end-i),2);      % Error signa, transformed domain
    eDr{i} = zeros(len,1);          % Error signal, time domain
    delays(i) = 2^i-1;              % Level delay for synthesis   
    w{i} = zeros(L(end-i),2);       % Subband adaptive filter coefficient, initialize to zeros    
    
end 
w{i} = zeros(L(end-i),3);           % Last level has 2 columns, cD and cA
eD{i} = zeros(1,3);                 % Last level has 2 columns, cD and cA
U.tmp = zeros(len,1);
Y.tmp = zeros(len,1);

%pwr = w;
%beta = 1./L(2:end-1);

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
    y = [dn(n); y(1:end-1)];        % Desired response vector        

    % Analysis Bank
    U.tmp = u;
    Y.tmp = y;
   
    
    for i = 1:level
        if mod(n,2^i) == 0    
            U.Z = H'*U.tmp;
            U.cA{i} = [U.Z(1); U.cA{i}(1:end-1)];
            U.cD1{i} = [U.Z(2); U.cD1{i}(1:end-1)]; 
            U.cD2{i} = [U.Z(3); U.cD2{i}(1:end-1)]; 
            U.tmp = U.cA{i}(1:len);
                    
            Y.Z = H'*Y.tmp;
            Y.cA{i} = [Y.Z(1); Y.cA{i}(1:end-1)];
            Y.cD1{i} = [Y.Z(2); Y.cD1{i}(1:end-1)]; 
            Y.cD2{i} = [Y.Z(3); Y.cD2{i}(1:end-1)]; 
            Y.tmp = Y.cA{i}(1:len);
             
            if i == level               
                filt_input =  [U.cA{i}, U.cD1{i},U.cD2{i} ];
                
                eD{i} = Y.Z' - sum((filt_input).*w{i});                    
                
                if n >= AdaptStart(i)
%                     pwr{i} = beta(i)*pwr{i}+ (1-beta(i))*([U.cA{i},U.cD{i}].*[U.cA{i},U.cD{i}]);
%                     w{i} = w{i} + mu*[U.cA{i},U.cD{i}].*eD{i}./((sum(pwr{i})+alpha)); 
                    w{i} = w{i} + mu*filt_input.*eD{i}./(sum(filt_input.*filt_input)+alpha); 
                                                          
                end 
            else                
                eD{i} = [eD{i}(2:end,:); [Y.cD1{i}(1), Y.cD2{i}(1)]  - sum([U.cD1{i}, U.cD2{i}].*w{i})];                 
                
                if n >= AdaptStart(i)
%                     pwr{i} = beta(i)*pwr{i}+ (1-beta(i))*(U.cD{i}.*U.cD{i});
%                     w{i} = w{i} + (mu*eD{i}(end)/((sum(pwr{i}) + alpha)))*U.cD{i};
                    w{i} = w{i} + (mu*eD{i}(end,:)./(sum([U.cD1{i}, U.cD2{i}].*[U.cD1{i}, U.cD2{i}])+ alpha)).*[U.cD1{i}, U.cD2{i}];                    
                    
                end
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
                eDr{i} = F*[eDr{i+1}(1); eD{i}((end-(len-1)*delays(end-i)),:)' ] + eDr{i}; %% problema qui cambiare
                eDr{i+1} = [eDr{i+1}(2:end); 0];                
              
            end            
        end
    end   
    en(n) = eDr{i}(1);
    eDr{i} = [eDr{i}(2:end); 0];        
end

en = en(1:ITER);
S.coeffs = w;
end