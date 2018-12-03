function [en,S] = DFThammerstein_adapt_malik(un,dn,S, L, select)
% DFT Subband Adaptive Hammerstein Filter (WAF)                 
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal
% L                 length of periodization for the non-linear part 
% select            input 'odd' 'even' or 'both' for both even and odd fourier coeffs 

% common parameters for both even and odd part 

power_alpha = S.palpha;        % Initialize constant for update bin power 
order = S.order; 
M = S.length;                  % Length of adaptive filter
mu_unconst = S.step_unconst;   % Step size of unconstrained FDAF
AdaptStart = 2*S.AdaptStart;   % Adaptation starts after 2M-sample blocks are acquired

dn_block = zeros(S.length,1);  % Initialize M samples of desired signal to 0

% ODD part 

if strcmp(select,'odd') || strcmp(select,'both')

WEIGHT_odd = zeros(2*S.length,order);  % Initialize 2M coefficients to 0 
u_odd_new = zeros(S.length,order);     % Initialize M samples of new block to 0
u_odd_old = zeros(S.length,order);     % Initialize M samples of old block to 0
power_odd = zeros(2*S.length,order);   % Initialize average power of each bin to 0

end

% EVEN part 

if strcmp(select,'even') || strcmp(select,'both')

WEIGHT_even = zeros(2*S.length,order);  % Initialize 2M coefficients to 0 
u_even_new = zeros(S.length,order);     % Initialize M samples of new block to 0
u_even_old = zeros(S.length,order);     % Initialize M samples of old block to 0
power_even = zeros(2*S.length,order);   % Initialize average power of each bin to 0

end

ITER = length(un);             % Length of input sequence      
en = [];                       % Accumulation of error vector 
alpha =S.alpha; 

for n = 1:ITER
    
    dn_block = [dn_block(2:end); dn(n)];  % Start desired signal block
    
    for i = 1:order
        
        if strcmp(select,'odd') || strcmp(select,'both')
        
            tmp_odd = sin(pi*(i)*un(n)/L);% cos(pi*i*un(n)/1);
            u_odd_new(:,i) = [u_odd_new(2:end,i); tmp_odd];        
        
        end
        
        if strcmp(select,'even') || strcmp(select,'both')
        
            tmp_even = cos(pi*(i)*un(n)/L);
            u_even_new(:,i) = [u_even_new(2:end,i); tmp_even]; 
        
        end
        
    end
    
    if mod(n,M)==0   
        % If iteration == block length, 
        
        for i = 1:order
            
            if strcmp(select,'odd') || strcmp(select,'both')
            
             un_blocks_odd(:,i) = [u_odd_old(:,i); u_odd_new(:,i)];       
             u_odd_old(:,i) = u_odd_new(:,i);
             
             end
             
             if strcmp(select,'even') || strcmp(select,'both')
             
             un_blocks_even(:,i) = [u_even_old(:,i); u_even_new(:,i)];       
             u_even_old(:,i) = u_even_new(:,i);
             
             end
        
        end
        
        if n >= AdaptStart                % Frequency-domain adaptive filtering
           
            format long

            % Set up a vector to extract the first M elements
 
            %window_save_first_M = [ones(1,M), zeros(1,M)]';  

            % Transform the reference signal into frequency domain
            
            yn = 0; 
            
            for i = 1:order
                
                if strcmp(select,'odd') || strcmp(select,'both')    
              
                 UN_odd(:,i) = fft(un_blocks_odd(:,i));
                 YN_odd(:,i) = UN_odd(:,i).*WEIGHT_odd(:,i);    % Multiplication of input and coeff. vectors
            
                end
            
                 if strcmp(select,'even') || strcmp(select,'both')  
            
                 UN_even(:,i) = fft(un_blocks_even(:,i));
                 YN_even(:,i) = UN_even(:,i).*WEIGHT_even(:,i);
            
                 end
                               
                                                        
            % Compute the estimate of desired signal
 
                 if  strcmp(select,'both')  
            
                     temp = real(ifft(YN_odd(:,i))) + real(ifft(YN_even(:,i))); % Real part of IFFT output
            
                 elseif strcmp(select,'odd') 
                 
                     temp = real(ifft(YN_odd(:,i))); 
                 
                 elseif   strcmp(select,'even') 
                 
                    temp = real(ifft(YN_even(:,i))); 
                 
                end
                 
                 
            
                yn = temp(M+1:2*M) + yn;             % Extracted the last M elements of IFFT block
            
             % Update the signal power estimates
            
            end
            
            % Compute the error signal in time domain

           error = dn_block-yn;              % Error signal block
           EN = fft([zeros(1,M),error']');   % Insert zero block to form 2M block before FFT

           % Compute the gradient and weight update in frequency domain 
           for i = 1:order
               
               if strcmp(select,'odd') || strcmp(select,'both')  
               
                   power_odd(:,i) = (power_alpha.*power_odd(:,i)) + (1 - power_alpha).*(abs(UN_odd(:,i)).^2); 
                   gradient_odd(:,i) = conj(UN_odd(:,i)).* EN;
                   WEIGHT_odd(:,i) = WEIGHT_odd(:,i) + mu_unconst.*(gradient_odd(:,i))./(power_odd(:,i)+ alpha);
               
               end
               
               
               if strcmp(select,'even') || strcmp(select,'both')  
                   
                   power_even(:,i) = (power_alpha.*power_even(:,i)) + (1 - power_alpha).*(abs(UN_even(:,i)).^2); 
                   gradient_even(:,i) = conj(UN_even(:,i)).* EN;
                   WEIGHT_even(:,i) = WEIGHT_even(:,i) + mu_unconst.*(gradient_even(:,i))./(power_even(:,i)+ alpha);
        
              
               end
             
               
           end
           
           en = [en; error];             % Update error block 
           
        end
        
            
        
                     
            
        end
end

if strcmp(select,'odd') || strcmp(select,'both')  
   S.weight_odd = WEIGHT_odd;
end

if strcmp(select,'even') || strcmp(select,'both')  
   S.weight_even = WEIGHT_even; 
end

en = en'; 

end


