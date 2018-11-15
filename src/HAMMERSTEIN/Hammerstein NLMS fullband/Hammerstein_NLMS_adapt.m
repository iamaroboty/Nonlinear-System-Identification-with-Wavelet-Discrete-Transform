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
p = zeros(order, 1); % non linearity coeffs vector 
p(1) = 1; 
X = zeros(M, order); 

ITER = length(un);              % Length of input sequence
en = zeros(1,ITER);             % Initialize error sequence to zero

	
for n = 1:ITER
    
    for i = 1:order % build the coeff vector (taylor expansion)        
         xp = [xp(2:end); un(n)^i];                   
    end 
    
    X = cat(1, xp', X(1:end-1,:));    
    
    
    % plot 
    
%     N = size(X,1);
%     n = (0:N-1)';
% 
% 
% 
% 
%     stem3(n,1:size(X,2),X','filled')
%     ax = gca;
%     ax.YTick = 1:3;
%     view(37.5,30)
% 
% [cr,lgs] = xcorr(X,'coeff');
%     
%     
% for row = 1:size(X,2)
%     for col = 1:size(X,2)
%         nm = size(X,2)*(row-1)+col;
%         subplot(size(X,2),size(X,2),nm)
%         stem(lgs,cr(:,nm),'.')
%         title(sprintf('c_{%d%d}',row,col))
%         ylim([0 1])
%     end
% end
    
      
    en(n) = dn(n) - w'*X*p;    
    
    if n >= AdaptStart       
        update_p = (mu(1)*en(n)*X'*w)/((X'*w)'*(X'*w) + alpha);
        update_w = (mu(2)*en(n)*X*p)/((X*p)'*(X*p) + alpha);
        p = (1-mu(1)*leak(1))*p + update_p;     
        w = (1-mu(2)*leak(2))*w + update_w;  

        S.iter = S.iter + 1;
    end
end

S.coeffs{1} = w;
S.coeffs{2} = p; 
end


    


