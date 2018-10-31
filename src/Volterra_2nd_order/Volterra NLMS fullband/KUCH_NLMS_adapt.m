function [en,S] = KUCH_NLMS_adapt(un,dn,S)
% Fullband Volterra filtering adapt, 2nd order              
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in Volterra_NLMS_init.m
% en                History of error signal

len_h1 = size(S.coeffs{1},1);                     % Adaptive filters length
Nc = size(S.coeffs{2},1); 
Nw = size(S.coeffs{3},1); 

mu1 = S.mu1;
muc = S.muc;
muw = S.muw;
AdaptStart = S.AdaptStart;        % Transient
alpha = S.alpha;                  % Small constant (1e-6)
%leak = S.leakage;                 % Leaky factor 

w1 = S.coeffs{1};              % Adaptive filter, linear contribute
c = S.coeffs{2};
w2 = S.coeffs{3};              % Adaptive filter, nonlinear contribute

x = zeros(len_h1,1);            % Delay line for first order input
%u = zeros(Nw, 1);  % Delay line for second order input

ITER = length(un);              % Length of input sequence
en = zeros(1,ITER);             % Initialize error sequence to zero


X = zeros(Nw, Nc);

	
for n = 1:ITER
    x = [un(n); x(1:end-1)];                    % First order input
    
    %u = [(c'*u1(1:S.coeffs{2}))^2; u(1:end-1)]; 
       
    X = cat(1, x(1:Nc)', X(1:end-1,:));
    
    v = X*c;
    u = v.^2;
    
    en(n) = dn(n) - w1'*x - w2'*u;
    
    
    
    if n >= AdaptStart
        w1 = w1 + (mu1*en(n)/(x'*x + alpha))*x; 
        
        norm_grade = (2*X'*diag(v)*w2)'*(2*X'*diag(v)*w2);
        c = c + muc*(en(n)/(norm_grade+ alpha))*(2*X'*diag(v)*w2);
        
        w2 = w2 + (muw*en(n)/(u'*u + alpha))*u; 
        S.iter = S.iter + 1;
    end
end

S.coeffs = [w1; c; w2];                     % Coefficient values at final iteration                 
end


    


