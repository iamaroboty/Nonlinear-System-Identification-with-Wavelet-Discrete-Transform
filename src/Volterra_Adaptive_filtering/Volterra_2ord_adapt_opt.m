function [en,S] = Volterra_2ord_adapt_opt(un,dn,S,C)
% Wavelet-Decomposition Subband Adaptive Filter (WAF)                 
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal

K = S.kernel_length;              % Kernel dimensions (1st and second order)
M = S.length;                     % filter lengths (array) 
mu = S.step;                      % Step Size here is an array 

AdaptStart = S.AdaptStart;        % Transient
alpha = S.alpha;                  % Small constant (1e-6)

% filter bank parameters
H = S.analysis;                   % Analysis filter bank
level = S.levels;                 % Wavelet Levels

H = create_multilevel_bank(H, level);
F = flip(H);

S.analysis_bank = H;
S.synthesis_bank = F;
[len, ~] = size(H);               % Wavelet filter length

if nargin == 3
    C = K(2);                 % Number of nonlinear channel
end

taplen = max(C,len);

U1_tot = zeros(M(1),2^level); % output of analysis filters for 1st order signal 

for i = 1:C
%     u2{i} = zeros(len, 1);      % Nonlinear input delay line
%     A2{i} = zeros(len,2^level);
    w2{i} = zeros(K(2)-i+1,1);
    
    U2_tot{i} = zeros(K(2)-i+1,2^level); % output of analysis filters for 2nd order signal 
end

u2 = zeros(len,C);
A2 = zeros(len,2^level,C);

d = zeros(len,1);               
u = zeros(taplen,1);
A1 = zeros(len,2^level);  

z = zeros(len,1);

w1 = zeros(M(1),1);


ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


for n = 1:ITER    
  
    d = [dn(n); d(1:end-1)];                        % Update tapped-delay line of d(n)
    u = [un(n); u(1:end-1)];                        % Update tapped-delay line of u(n)

    A1 = [u(1:len), A1(:,1:end-1)];                 % Update buffer linear signal
    
    % Building delay line for quadratic adaptive filter
    for i = 1:C
        u2(:,i) = [un(n)*u(i); u2(1:end-1,i)];
        A2(:,:,i) = [u2(1:len,i), A2(:,1:end-1,i)];       % use 3d matrix
    end
  
       
    if (mod(n,2^level)==0)                               % Tap-weight adaptation at decimated rate
        U1 = (H'*A1)';                              % Partitioning u(n) 
        U1_2 = U1_tot(1:end-2^level,:);
        U1_tot = [U1', U1_2']';                           % Subband data matrix
        
        for i = 1:C
            U2{i} = (H'*A2(:,:,i))';                              % Partitioning u(n) 
            U2_2{i} = U2_tot{i}(1:end-2^level,:);
            U2_tot{i} = [U2{i}', U2_2{i}']';                           % Subband data matrix
        end
        
        dD = H'*d;                                 % Partitioning d(n) 
        
        e2 = 0;
        norm = 0;
        for i = 1:C
            
            e2 = e2 +U2_tot{i}(1:size(w2{i},1),:)'*w2{i};
            
            norm = norm + sum(U2_tot{i}(1:size(w2{i},1),:).*U2_tot{i}(1:size(w2{i},1),:));
        end
        
        eD = dD - U1_tot'*w1 - e2;                         % Error estimation
        
        if n >= AdaptStart             
            w1 = w1 + U1_tot*(eD./(sum(U1_tot.*U1_tot)+alpha)')*mu(1); % Tap-weight adaptation          
            
            for i= 1:C
                w2{i} = w2{i} + U2_tot{i}(1:size(w2{i},1),:)*(eD./(norm +alpha)')*mu(2);
            end
           
        end
        z = F*eD + z;                                       
%         en(n-2^level+1:n) = z(1:2^level); 
%         z = [z(2^level+1:end); zeros(2^level,1)]; 
    end
    en(n) = z(1);
    z = [z(2:end); 0];                         
end


S.coeffs{1}=w1; 
S.coeffs{2}=w2;
end


    


