function [en,S] = PHammerstein_NLMS_adapt(un,dn,S)
% Fullband Volterra filtering adapt, 2nd order              
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in Volterra_NLMS_init.m
% en                History of error signal

n_filters = size(S.filters_lengths,2); 
M = S.filters_lengths;              % kernel memory lengths 
mu = S.step;                      % Step Size here is an array 
AdaptStart = S.AdaptStart;        % Transient
alpha = S.alpha;                  % Small constant (1e-6)
leak = S.leaks;                 % Leaky factor 


for i = 1:n_filters
    u{i} = zeros(M(i), 1);      
    w{i} = zeros(M(i),1);
end


ITER = length(un);              % Length of input sequence
en = zeros(1,ITER);             % Initialize error sequence to zero

	
for n = 1:ITER
    
    sum = 0 ; 
    norm = 0; 
    for i = 1:n_filters
         u{i} = [un(n)^i; u{i}(1:end-1)];                    % First order input
         sum = sum + w{i}'*u{i};
         norm = norm + u{i}'*u{i};       
    end 
    en(n) = dn(n) - sum;    
    
    if n >= AdaptStart       
        for i = 1:n_filters
            w{i} = (1-mu(i)*leak(i))*w{i} + (mu(i)*en(n)/(norm + alpha))*u{i}; 
        end
        S.iter = S.iter + 1;
    end
end

S.coeffs = w;           
end


    


