function [out] = Hammerstein_NLMS_test(un,S)            
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
w = S.coeffs{1};
p = S.coeffs{2};

X = zeros(M, order); 

ITER = length(un);              % Length of input sequence
out = zeros(1,ITER);             % Initialize error sequence to zero

	
for n = 1:ITER
    for i = 1:order % build the coeff vector (taylor expansion)        
         xp = [xp(2:end); un(n)^i];                   
    end 
    
    X = cat(1, xp', X(1:end-1,:));   
      
    out(n) = w'*X*p;    
end

end


    


