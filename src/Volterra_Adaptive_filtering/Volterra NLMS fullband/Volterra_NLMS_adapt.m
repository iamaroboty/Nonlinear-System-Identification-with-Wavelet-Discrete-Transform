function [en,S] = Volterra_NLMS_adapt(un,dn,S)
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

M = S.length;                     % Adaptive filters length
mu = S.step;                      % Step Size here is an array 
AdaptStart = S.AdaptStart;        % Transient
alpha = S.alpha;                  % Small constant (1e-6)
leak = S.leakage;                 % Leaky factor 

w1 = zeros(M(1),1);              % Adaptive filter, linear contribute
w2 = zeros(M(2),1);              % Adaptive filter, nonlinear contribute

u1 = zeros(M(1),1);              % Delay line for first order input
utmp = zeros(K(2)^2-K(2)+1,1);  % Temporal Delay line for second order input - maybe not needed
u2 = zeros(M(2),1);             % Delay line for second order input

index = [];                     % Indexing for correct second order input

for i=0:K(2)-1
    index = cat(2,index, K(2)*i+1:K(2)*(i+1)-i);
end

ITER = length(un);              % Length of input sequence
en = zeros(1,ITER);             % Initialize error sequence to zero

	
for n = 1:ITER
    u1 = [un(n); u1(1:end-1)];                    % First order input
    utmp = [un(n)*u1(1:K(2)); utmp(1:end-K(2))];    % Temporal delay line
    u2 = utmp(index);                           % Second order input
    
    en(n) = dn(n) - w1'*u1 - w2'*u2;
    
    if n >= AdaptStart
        w1 = (1-mu(1)*leak)*w1 + (mu(1)*en(n)/(u1'*u1 + alpha))*u1; 
        w2 = (1-mu(2)*leak)*w2 + (mu(2)*en(n)/(u2'*u2 + alpha))*u2; 
        S.iter = S.iter + 1;
    end
end

S.coeffs = {w1, w2};                     % Coefficient values at final iteration                 
end


    


