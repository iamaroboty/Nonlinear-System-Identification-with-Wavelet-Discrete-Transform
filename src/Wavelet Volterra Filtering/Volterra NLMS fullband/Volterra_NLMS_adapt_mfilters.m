function [en,S] = Volterra_NLMS_adapt_2(un,dn,S,C)
% Fullband Volterra filtering adapt, 2nd order              
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in Volterra_NLMS_init.m
% en                History of error signal

K = S.kernel_length;              % kernel memory lengths 

if K(2) > K(1)
    error('2nd order kernel memory must be shorter than 1st one\n');
end

if nargin == 3
    C = K(2) ;
end

M = S.length;                     % Adaptive filters length
mu = S.step;                      % Step Size here is an array 
AdaptStart = S.AdaptStart;        % Transient
alpha = S.alpha;                  % Small constant (1e-6)
leak = S.leakage;                 % Leaky factor 

w1 = zeros(M(1),1);              % Adaptive filter, linear contribute

% no indexing but separate filtering
taplen = C;
u = zeros(M(1),1);              % Delay line for first order input

for i = 1:C
    u2{i} = zeros(K(2)-i+1, 1);      
    w2{i} = zeros(K(2)-i+1,1);
end


ITER = length(un);              % Length of input sequence
en = zeros(1,ITER);             % Initialize error sequence to zero

	
for n = 1:ITER
    u = [un(n); u(1:end-1)];                    % First order input
    
    sum = 0; 
    norm = 0; 
    for i = 1:C
        u2{i} = [un(n)*u(i); u2{i}(1:end-1)];    % second order input 
        sum = sum + w2{i}'*u2{i};
        norm = norm + u2{i}'*u2{i};
        
    end  
    norm = norm + u'*u;
    
    en(n) = dn(n) - w1'*u(1:M(1)) - sum;
    
    if n >= AdaptStart
        w1 = (1-mu(1)*leak)*w1 + (mu(1)*en(n)/(norm + alpha))*u; 
    for i = 1:C
        w2{i} = (1-mu(2)*leak)*w2{i} + (mu(2)*en(n)/(norm + alpha))*u2{i}; 
    end
        S.iter = S.iter + 1;
    end
end

S.coeffs{1} = w1;
S.coeffs{2} = w2;                     % Coefficient values at final iteration                 
end


    


