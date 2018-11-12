function [en,S] = Hammerstein_NLMS_adapt(un,dn,S)            
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters 
% en                History of error signal

order = S.order; 
M = S.filters_lengths;              % kernel memory lengths 
mu = S.step;                      % Step Size here is an array 
AdaptStart = S.AdaptStart;        % Transient
alpha = S.alpha;                  % Small constant 
leak = S.leaks;                 % Leaky factor 



xp = zeros(order, 1);  
w = zeros(M,1);
w(1) = sqrt(order);
p = zeros(order, 1); % non linearity coeffs vector 
p(1) = sqrt(order); 
X = zeros(M, order); 

ITER = length(un);              % Length of input sequence
en = zeros(1,ITER);             % Initialize error sequence to zero

	
for n = 1:ITER
    
    for i = 1:order % build the coeff vector (taylor expansion)        
         xp = [xp(2:end); un(n)^i];                   
    end 
    
    X = cat(1, xp', X(1:end-1,:));      
    
      
    en(n) = dn(n) - w'*X*p;    
    
    if n >= AdaptStart        
        p = (1-mu(1)*leak(1))*p + (mu(1)*en(n)*X'*w)/((X'*w)'*(X'*w) + alpha);        
        w = (1-mu(2)*leak(2))*w + (mu(2)*en(n)*X*p)/((X*p)'*(X*p) + alpha);  

        S.iter = S.iter + 1;
    end
end

S.coeffs{1} = w;
S.coeffs{2} = p; 
end


    


