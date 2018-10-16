function [en,S] = Adapt_2layer_af(un,dn,S)
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


M = S.length;                     % Unknown system length (Equivalent adpative filter lenght)

F = S.analysis;                   % Analysis filter bank
H = S.synthesis;                  % Synthesis filter bank 


% petraglia aliasing free structure adaptation

% filters for the aliasing free bank 

%upsampled filters


% %check for two layers
check_H_extd = cat(2, conv(H(:,1), upsample(H(:,1),2)), conv(H(:,1), upsample(H(:,2),2)), conv(H(:,2), upsample(H(:,1),2)), conv(H(:,2), upsample(H(:,2),2)) ); 
% H0, H1, H2, H3, H4
% Hi = upsample(H,2);
% Hi = [conv(Hi(:,1),H(:,1)), conv(Hi(:,2),H(:,1)), conv(Hi(:,1),H(:,2)), conv(Hi(:,2),H(:,2))];  
% if mod(length(Hi),2) ~= 0
%     Hi = Hi(1:end-1,:);
% end
% S.analysis = Hi;
% S.synthesis = flip(Hi);
% 

Hi = check_H_extd(1:end-1,:);


H_af = cat(2, conv(check_H_extd(:,1), check_H_extd(:,1)), conv(check_H_extd(:,1), check_H_extd(:,2)), ...
                conv(check_H_extd(:,2), check_H_extd(:,2)),  conv(check_H_extd(:,2), check_H_extd(:,3)), ...
                  conv(check_H_extd(:,3), check_H_extd(:,3)),  conv(check_H_extd(:,3), check_H_extd(:,4)), ...
                   conv(check_H_extd(:,4), check_H_extd(:,4))   ); 
               
        

F = flip(Hi); 


% analysis and synthesis are used in reverse to obtain in U.Z a column
% vector with cD in the first position

[len_af, ~] = size(H_af);               % Wavelet filter length
[len, ~] = size(Hi); 

level = S.levels;                 % Wavelet Levels
L = S.L;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

% Init Arrays
% everything is brought to the first level
    
U_c = zeros(L(end-level),2^(level+1)-1);  
Y_c = zeros(L(end-level),2^(level)); 

eDr = zeros(len,1);          % Error signal, time domain
delay = 1;                    % Level delay for synthesis
           
w = zeros(L(end-level),2^level);           % Last level has 2 columns, cD and cA



eD = zeros(1,2^level);              % Last level has 2 columns, cD and cA

pwr = w;
beta = 1./L(2:end-1);

u = zeros(len_af,1);                 % Tapped-delay line of input signal (Analysis FB)  
y = zeros(len,1); 

ytap = zeros(len_af-len,1); 

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    

    y = [ytap(end); y(1:end-1)];        % Desired response vector
    
    ytap = [dn(n); ytap(1:end-1)]; 
    
    
    U.tmp = u;
    Y.tmp = y; 
    
        if mod(n,2^level) == 0
                       
            U.Z = H_af'*U.tmp; % column [cD ; cA] 
            Y.Z = Hi'*Y.tmp; 
                     
            [rows, cols] = size(U.Z);
            
            indx = 1;
            
            for col=1:cols
                for row=1:rows 
                    
                    U_c(:,indx) = cat(1,U.Z(row,col), U_c(1:end-1, indx)); %CD||CA
        
                indx=indx+1;
                end  
            end
                          
%             direct = zeros(1 ,2^level); 
%             
%             indx = 1; 
%             
%             % direct nodes 
%             for j=1:2:size(U_c,2)
%             direct(:,indx) = sum(U_c(:,j).*w(:,indx));
%             indx = indx +1; 
%             end
%             
%             cross = zeros(1 ,2^(level+1)-2);
%             
%             indx1 = 1; 
%             indx2 = 2; 
%             
%             %cross nodes 
%             for j=2:2:size(U_c,2)
%             cross(:,indx1) = sum(U_c(:,j).*w(:,indx2));
%             indx1 = indx1+1; 
%             indx2 = indx2 -1; 
%             cross(:,indx1) = sum(U_c(:,j).*w(:,indx2));
%             indx1 = indx1+1; 
%             indx2 = indx2 +2; 
%             end
%             
%             % sum nodes 
%             tmp = zeros(1 ,2^level); 
%             
%              
%             tmp(:,1) = cross(:,1);
%             indx = 2;
%             
%             for j=2:2:size(cross,2)-1
%                tmp(:,indx) = cross(:,j) + cross(:, j+1);
%                indx = indx+1; 
%               
%                 
%             end
%             
%             tmp(:,end) = cross(:,end);

            directH0H0 = sum(U_c(:,1).*w(:,1)); 
            directH1H1 = sum(U_c(:,3).*w(:,2)); 
            directH2H2 = sum(U_c(:,5).*w(:,3)); 
            directH3H3 = sum(U_c(:,7).*w(:,4)); 
            
            crossH1H0G1 = sum(U_c(:,2).*w(:,2));
            crossH1H0G0 = sum(U_c(:,2).*w(:,1));
            crossH2H1G2 = sum(U_c(:,4).*w(:,3));
            crossH2H1G1 = sum(U_c(:,4).*w(:,2));
            crossH3H2G3 = sum(U_c(:,6).*w(:,4));
            crossH3H2G2 = sum(U_c(:,6).*w(:,2));
            
            summed = [directH0H0+crossH1H0G1; directH1H1+crossH1H0G0+crossH2H1G2; directH2H2+crossH2H1G1+crossH3H2G3; ...
                   directH3H3+crossH3H2G2];
            
            eD = Y.Z - (summed) ;
            
            
            if n >= AdaptStart

                    
            w = w + mu.*(U_c(:,1:2:end).*eD'./(sum(U_c(:,1:2:end).*U_c(:,1:2:end))+alpha)); 
                                  
            end       
            
            % Synthesis 
            eDr = F*eD ;
           
            S.iter{1} = S.iter{1} + 1;  
            
            en(n-2^level+1:n) = eDr(1:2^level);
            
        end
            
      
end
en = en(1:ITER);
S.coeffs = w;
end


