function [en,S] = DFThammerstein_adapt_malik(un,dn,S)
% DFT Subband Adaptive Hammerstein Filter (WAF)                 
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal


% mu = S.step;                      % Step Size
% AdaptStart = S.AdaptStart;        % Transient
% alpha = S.alpha;                  % Small constant (1e-6)
% 
% 
% %H = H(:,1:N_subbands/2+1); % Analysis filter (only the first N/2+1 are taken)
% 
% order = S.order; 
% 
% W = S.coeffs{1}; 
% % p(1) = 1; 
% 
% %C = zeros(size(w,1),size(w,2), order);                 
% %a = zeros(len,1); d = zeros(len,1); A = zeros(len, size(w,2));  
% % z = zeros(len,1);
% f = zeros(size(W,1),order); d = zeros(size(W,1),1);
% 
% F = S.dftmat; 
% 
% G = F*F';
% 
% ITER = length(un);
% en = zeros(1,ITER);    % Initialize error sequence to zero
% 
% for n = 1:ITER
%     
%     d = [ dn(n); d(1:end-1)]; 
%     
%     for j= 1:order
%         
%        %tmp_d = sin(pi*j*dn(n)/100)+ cos(pi*j*dn(n)/100);
%        tmp_u = sin(pi*j*un(n)/100)+ cos(pi*j*dn(n)/100);   
%                            % Update tapped-delay line of d(n)
%                            
%        f(:,j) = [ tmp_u; f(1:end-1,j)]; % Update tapped-delay line of u(n)
%        
%        X(:,:,j) = diag(S.dftmat*f(:,j)); 
%     
%     end
%     
%     if mod(n, size(W,1))
%        
%         
%         
%         
%     end
%       
%     Y = diag(S.dftmat*d); 
%     
%     E = Y - G*X*W; 
%     
%     if n > S.AdaptStart
%         
%        
%          
%     end
%     %if (mod(n-1,decfac)==0)                         % Tap-weight adaptation at decimated rate

%       
% end
% 
% 
% end

order = S.order; 
M = S.length;                  % Length of adaptive filter
mu_unconst = S.step_unconst;   % Step size of unconstrained FDAF
AdaptStart = 2*S.AdaptStart;   % Adaptation starts after 2M-sample blocks are acquired
WEIGHT = zeros(S.length,order);  % Initialize 2M coefficients to 0 
dn_block = zeros(S.length,1);  % Initialize M samples of desired signal to 0
u_new = zeros(S.length,order);     % Initialize M samples of new block to 0
power_alpha = S.palpha;        % Initialize constant for update bin power 
power = zeros(S.length,order);   % Initialize average power of each bin to 0
ITER = length(un);             % Length of input sequence      
en = [];                       % Accumulation of error vector 
 

F = dftmtx(M); 
Q = vertcat(zeros(M,M),eye(M));  
G = F*(Q')*Q*F'; 

for n = 1:ITER
    
    dn_block = [dn_block(2:end); dn(n)];  % Start desired signal block
    
    for i = 1:order
        tmp = sin(pi*i*un(n)/100) + cos(pi*i*un(n)/100);
        u_new(:,i) = [u_new(2:end,i); tmp];        % Start input signal blocks    
        
    end
    
    if mod(n,M)==0   
        % If iteration == block length, 
        
        for i = 1:order
            
             un_blocks(:,i) = [ u_new(:,i)];       
             %u_old(:,i) = u_new(:,i);
        
        end
        
        if n >= AdaptStart                % Frequency-domain adaptive filtering
           
            format long

            % Set up a vector to extract the first M elements
 
            %window_save_first_M = [ones(1,M), zeros(1,M)]';  

            % Transform the reference signal into frequency domain
            
            %yn = 0; 
            
            Yn = F*dn_block; 
            
            for i = 1:order
            X(:,:,i) = diag(F*(un_blocks(:,i)));
            C(:,:,i) = G*X(:,:,i); 
                                  % FFT[old block; new block]
                                  % Old block contains M previous input samples (u_old)
                                  % New block contains M new input samples (u_new)
                                                        
            % Compute the estimate of desired signal
 
            %YN(:,i) = UN(:,i).*WEIGHT(:,i);       % Multiplication of input and coeff. vectors
            %temp = real(ifft(YN(:,i)));       % Real part of IFFT output
            %yn = temp(M+1:2*M) + yn;             % Extracted the last M elements of IFFT block
            EN = Yn - C(:,:,i)*WEIGHT(:,i); 
            
            power(:,i) = (power_alpha.*power(:,i)) + (1 - power_alpha).*diag((ctranspose(X(:,:,i))*X(:,:,i))); 
            
            gradient(:,i) = diag(ctranspose(X(:,:,i))).*EN;
            WEIGHT(i) = WEIGHT(i) + mu_unconst.*(gradient(i));
            
            end
            
            error = ifft(EN); 

            % Compute the error signal

           %error = dn_block-yn;              % Error signal block
           %EN = fft([zeros(1,M),error']');   % Insert zero block to form 2M block before FFT


           % Update the signal power estimates

           

           % Compute the gradient and weight update in frequency domain 

          
           
           
          
        
            
        en = [en; error];             % Update error block 
                     
            
        end
    end
end

S.weight = WEIGHT;
en = en'; 

end


