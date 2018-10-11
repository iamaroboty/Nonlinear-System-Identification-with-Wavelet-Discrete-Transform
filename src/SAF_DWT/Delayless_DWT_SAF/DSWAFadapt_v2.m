function [en,S] = DSWAFadapt_v2(un,dn,S)
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
[len, ~] = size(H);               % Wavelet filter length
level = S.levels;                 % Wavelet Levels
L = S.L;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]
UpdateRate = S.UpdateRate;
Hadam = hadamard(2^level);

% Init Arrays
for i= 1:level
       U.c{i} = zeros(L(end-i),2^(i));           
       E.c{i} = zeros(L(end-i),2^(i));      % Error signa, transformed domain
%        delays(i) = 2^i-1;                    % Level delay for synthesis    
       %w{i} = zeros(L(end-i),1);       % Subband adaptive filter coefficient, initialize to zeros    
end 
w = zeros(L(end-i),2^i);           % Last level has 2 columns, cD and cA
W = zeros(1,M);  
U.tmp = zeros(len,1);
E.tmp = zeros(len,1);
U.Z = zeros(2,1); %!!
E.Z = zeros(2,1); %!!
% pwr = w;
% beta = 1./L(2:end-1);

x = zeros(M,1);                    % Fullband tap-input vector
u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  
e = zeros(len,1);                 % Tapped-delay line of error

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero

if isfield(S,'unknownsys')
    b = S.unknownsys;
    norm_b = norm(b);
    eml = zeros(1,ITER);
    ComputeEML = 1;
else
    ComputeEML = 0;
end


for n = 1:ITER    
    x = [un(n); x(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    u = [un(n); u(1:end-1)];        % Desired response vector  

    if ComputeEML == 1
        eml(n) = norm(b-W')/norm_b; % System error norm (normalized)
    end    
    
    en(n) = dn(n) - W*x;         % Error signal
    e = [en(n); e(1:end-1)];
    
    U.tmp = u;
    E.tmp = e;
    
    for i = 1:level
        if mod(n,2^i) == 0                       
            U.Z = H'*U.tmp;
            E.Z = H'*E.tmp;
            
            [rows, cols] = size(U.Z);            
            indx = 1;
            
            for col=1:cols
                for row=1:rows                     
                    U.c{i}(:,indx) = cat(1,U.Z(row,col), U.c{i}(1:end-1, indx));
                    E.c{i}(:,indx) = cat(1,E.Z(row,col), E.c{i}(1:end-1, indx));
              
                    indx=indx+1;
                end  
            end
            
            U.tmp = U.c{i}(1:len,:);   
            E.tmp = E.c{i}(1:len,:);
            
            if i == level                                                
                if n >= AdaptStart
%                     pwr{i} = beta(i)*pwr{i}+ (1-beta(i))*([U.cA{i},U.cD{i}].*[U.cA{i},U.cD{i}]);
%                     w{i} = w{i} + mu*[U.cA{i},U.cD{i}].*eD{i}./((sum(pwr{i})+alpha)); 
                      w = w + mu*U.c{i}.*E.c{i}(1,:)./(sum(U.c{i}.*U.c{i})+alpha);                                                 
                end                     
            end           
            S.iter{i} = S.iter{i} + 1;                
        end
    end 
    
    if (n >= AdaptStart+M) && (mod(n,UpdateRate)==0)          
%         W = (waverec([w{i}(:,1); w{i}(:,2)], L, F(:,2), F(:,1)))';
        if level == 1
            W = 1/2.* Hadam * [w(:,2), w(:,1)]';
            W = W(:)'; 
        elseif level == 2
            ww = [w(:,4), w(:,2), w(:,3), w(:,1)];
            W = 1/4 .* Hadam * ww';
            W = W(:)'; 
%             W = flip(W);
        end              
    end
         
end

S.FULLcoeffs = W;
S.SUBcoeffs = w;
if ComputeEML == 1
    S.eml = eml;
end


