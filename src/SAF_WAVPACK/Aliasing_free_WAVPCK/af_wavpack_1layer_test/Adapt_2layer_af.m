function [en,S] = Adapt_2layer_af(un, dn, S)
% Wavelet-Decomposition Subband Adaptive Filter (WAF)                 
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal

M = S.length;                     % Unknown system length (Equivalent adpative filter lenght)
mu = S.step;                      % Step Size
AdaptStart = S.AdaptStart;        % Transient
alpha = S.alpha;                  % Small constant (1e-6)

F = S.analysis;                   % Analysis filter bank
H = S.synthesis;                  % Synthesis filter bank
level = S.levels;                 % Wavelet Levels
L = S.L;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

Hi = zeros(2^(level-1)*size(H,1), 2^(level));

indx = 1;



for i = 1:size(H,2)
for j=1:level-1
    
   up{i,j} = upsample(H(:,i), 2^(j)); 
   
end
end


%outer product

% Hi equivalent one level filters

H_tmp = H; 
for i=1:size(up,2)

H_tmp = outer_conv(H_tmp, up(:,i));
end

Hi = H_tmp; 

if mod(length(Hi),2) ~= 0
    Hi = Hi(1:end-1,:);
end

tmp = Hi(:,3); 
Hi(:,3) = Hi(:,end);
Hi(:,end) = tmp; 

% petraglia's structure af filters

indx = 1; 
indx2 = 1 ; 

for i= 1:size(Hi,2)
                          
      H_af(:,indx) = conv(Hi(:,i), Hi(:,i));
      indx = indx +1; 
      
      if i+1 <= size(Hi,2)
      H_af(:,indx) = conv(Hi(:,i), Hi(:,i+1));
      indx = indx +1; 
      end         
    
end



F = flip(Hi); 

%% petraglia aliasing free structure adaptation



% analysis and synthesis are used in reverse to obtain in U.Z a column
% vector with cD in the first position

[len_af, ~] = size(H_af);         % Wavelet filter length
[len, ~] = size(Hi); 

level = S.levels;                 % Wavelet Levels
L = S.L;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

% Init Arrays

% everything is brought to the first level
    
U_c = zeros(L(end-level),2^(level+1)-1);            
eDr = zeros(len,1);          % Error signal, time domain
z = zeros(len,1);            % reconstructed error signal
           
w = zeros(L(end-level),2^level);  % 2^level filters

u = zeros(len_af,1);              % Tapped-delay line of input signal (Analysis FB)  
y = zeros(len,1);                 % Tapped-delay line of desired response (Analysis FB)

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero
eDvec = zeros(2^level, ITER);

ytap = zeros(len_af-len,1); 
yztap = zeros(2^level, floor((len_af+len)/8));

for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    y = [dn(n); y(1:end-1)]; 
   
    %ytap = [dn(n); ytap(1:end-1)]; % tapped delay line for desired signal 
    
    U.tmp = u;
    Y.tmp = y; 
    
        if (mod(n,2^level) == 0)
            
            U.Z = H_af'*U.tmp;
            Y.Z = Hi'*Y.tmp;
%             Y.Z = yztap(:,end); %Hi'*Y.tmp;
%             yztap = cat(2, Hi'*Y.tmp, yztap(:,1:end-1));
         
            [rows, cols] = size(U.Z);
            
            indx = 1;
            
            for col=1:cols
                for row=1:rows 
                    
                    U_c(:,indx) = cat(1,U.Z(row,col), U_c(1:end-1, indx)); %CD||CA
        
                indx=indx+1;
                end  
            end
            
            directH0H0 = sum(U_c(:,1).*w(:,1)); 
            directH1H1 = sum(U_c(:,3).*w(:,2)); 
            directH2H2 = sum(U_c(:,5).*w(:,3)); 
            directH3H3 = sum(U_c(:,7).*w(:,4)); 
            
            crossH1H0G1 = sum(U_c(:,2).*w(:,2));
            crossH1H0G0 = sum(U_c(:,2).*w(:,1));
            crossH2H1G2 = sum(U_c(:,4).*w(:,3));
            crossH2H1G1 = sum(U_c(:,4).*w(:,2));
            crossH3H2G3 = sum(U_c(:,6).*w(:,4));
            crossH3H2G2 = sum(U_c(:,6).*w(:,3));
            
            summed = [directH0H0+crossH1H0G1; directH1H1+crossH1H0G0+crossH2H1G2;...
                        directH2H2+crossH2H1G1+crossH3H2G3; ...
                        directH3H3+crossH3H2G2];
            
             eD =  Y.Z - (summed) ;  
             eDvec(:,n) =eD;
             
             if n >= AdaptStart
                 
             w = w + mu.*(U_c(:,1:2:end).*eD'./(sum(U_c(:,1:2:end).*U_c(:,1:2:end))+alpha));              
             
             end

            
           z = F*eD + z;                                       
           en(n-2^level+1:n) = z(1:2^level); 
           z = [z(2^level+1:end); zeros(2^level,1)]; 
            
        end
      
end
en = en(1:ITER);
S.coeffs = w;
end


