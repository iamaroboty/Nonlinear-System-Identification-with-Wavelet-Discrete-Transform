%% WAVELET PACKET DECOMPOSITION TEST FOR PERFECT RECONSTRUTION

addpath 'Common';             % Functions in Common folder
% clear all; close all

% Testing Signal

d = 256;        %Total signal length
t=0:0.001:10;
un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
un = un(1:d)';
Ovr = 1; 


%% wavpack parameters

mu = 0.1;                      % ignored here 
M = 256;                        % Length of unknown system response also ignored here
level = 2;                     % Levels of Wavelet decomposition
filters = 'db1';               % Set wavelet type


S = QMFInit(M,mu,level,filters);
%S = SWAFinit(M, mu, level, filters); 

M = S.length;                     % Unknown system length (Equivalent adpative filter lenght)

F = S.analysis;                   % Analysis filter bank
H = S.synthesis;                  % Synthesis filter bank 


%% petraglia aliasing free structure adaptation

% filters for the aliasing free bank 

%upsampled filters


% %check for two layers
check_H_extd = cat(2, conv(H(:,1), upsample(H(:,1),2)), conv(H(:,1), upsample(H(:,2),2)), conv(H(:,2), upsample(H(:,1),2)), conv(H(:,2), upsample(H(:,2),2)) ); 
% H0, H1, H2, H3, H4
Hi = upsample(H,2);
Hi = [conv(Hi(:,1),H(:,1)), conv(Hi(:,2),H(:,1)), conv(Hi(:,1),H(:,2)), conv(Hi(:,2),H(:,2))];  
if mod(length(Hi),2) ~= 0
    Hi = Hi(1:end-1,:);
end
S.analysis = Hi;
S.synthesis = flip(Hi);

check_H_extd = check_H_extd(1:end-1,:);


check_H_af = cat(2, conv(check_H_extd(:,1), check_H_extd(:,1)), conv(check_H_extd(:,1), check_H_extd(:,2)), ...
                conv(check_H_extd(:,2), check_H_extd(:,2)),  conv(check_H_extd(:,2), check_H_extd(:,3)), ...
                  conv(check_H_extd(:,3), check_H_extd(:,3)),  conv(check_H_extd(:,3), check_H_extd(:,4)), ...
                   conv(check_H_extd(:,4), check_H_extd(:,4))   ); 

H_af = check_H_af;
H = check_H_extd; 
%F = cat(2, conv(F(:,1), upsample(F(:,1),2)), conv(F(:,1), upsample(F(:,2),2)), conv(F(:,2), upsample(F(:,1),2)), conv(F(:,2), upsample(F(:,2),2)) ); 
F = flip(H);   

%check one layer
% H_af = cat(2, conv(H(:,1), H(:,1)), conv(H(:,1), H(:,2)), conv(H(:,2), H(:,2))); 


% aliasing free petraglia filters

% H_af = zeros(2^(level)*size(H,1)-1, 2^(level+1)-1);
% col = 1; 
% already_computed = 0; 
% 


% for i =1:size(H,2) 
%     
% for j=1:size(H,2)
%     
%     
%     temp = H(:,i); 
%     
%     if already_computed == 1
%          already_computed = 0;
%         continue; 
%     end
%     
%     for indx=1:level
%     
%     temp = conv(temp, upsample(H(:,j), 2^(indx-1)));
%  
%         
%     end
%     
%     if j-i==1
%        already_computed =1;  
%     end
%     
%     H_af(:,col) = temp;
%     col = col +1;   
%     
% end
% 
% end
% 
% % equivalent filters in one level: 
% 
% 
% H_lvl = zeros(2^(level-1)*size(H,1), 2^(level));
% 
% for i=1:size(H(:,2))
%     
%     tmp=H(:,i); 
% 
% for j=2:level
%     
%     tmp = conv(tmp, upsample(H(:,i), 2^(j-1))); 
%     
% 
% end
% 
% H_lvl(:,i) = tmp; 
% 
% 
% end




% analysis and synthesis are used in reverse to obtain in U.Z a column
% vector with cD in the first position


[len_af, ~] = size(H_af);               % Wavelet filter length
[len, ~] = size(H); 

level = S.levels;                 % Wavelet Levels
L = S.L.*Ovr;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

% Init Arrays

% everything is brought to the first level
    
U_c = zeros(L(end-level),2^level);            
eDr = zeros(len,1);          % Error signal, time domain
delay = 1;                    % Level delay for synthesis
           
w = zeros(L(end-level),2^level);           % Last level has 2 columns, cD and cA

w(1,:) = 1;                   % set filters to kronecker delta

eD = zeros(1,2^level);              % Last level has 2 columns, cD and cA

pwr = w;
beta = 1./L(2:end-1);

u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'

    % Analysis Bank
    U.tmp = u;
    
        if mod(n,2^level) == 0
            
            
            U.Z = H'*U.tmp; % column [cD ; cA] 
            
         
            [rows, cols] = size(U.Z);
            
            indx = 1;
            
            for col=1:cols
                for row=1:rows 
                    
                    U_c(:,indx) = cat(1,U.Z(row,col), U_c(1:end-1, indx)); %CD||CA
        
                indx=indx+1;
                end  
            end
            
            %U.tmp = U_c(1:len,:);    
            
            direct = zeros(1 ,2^level); 
            
            indx = 1; 
            
            % direct nodes 
            for j=1:1:size(U_c,2)
            direct(:,indx) = sum(U_c(:,j).*w(:,indx));
            indx = indx +1; 
            end
            
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
%               
%                 
%             end
%             
%             tmp(:,end) = cross(:,end);
            
            
            
                
                eD = direct' ; 
                
                
                
            % Synthesis 
            eDr = F*eD ;
                     
           
            S.iter{1} = S.iter{1} + 1;  
            
            en(n-2^level+1:n) = eDr(1:2^level);
            
        end
        
        
        
      
end
en = en(1:ITER);


%% check for perfect reconstruction

tot_delay = 1;%(2^level - 1)*(len-1) +1 ;

stem(en(tot_delay:end));
hold on;
stem(un); 

%% Full packet 2 layer : TESTING
% % Init Arrays
% for i= 1:level
%     U.cD{i} = zeros(L(end-i),1);    
%     U.cA{i} = zeros(L(end-i),1);    
%     U.cD2{i} = zeros(L(end-i),1);    
%     U.cA2{i} = zeros(L(end-i),1); 
%     Y.cD{i} = zeros(L(end-i),1);
%     Y.cA{i} = zeros(L(end-i),1);
%     Y.cD2{i} = zeros(L(end-i),1);
%     Y.cA2{i} = zeros(L(end-i),1);    
%     eD{i} = zeros(L(end-i),1);      % Error signa, transformed domain
%     eDr{i} = zeros(len,1);          % Error signal, time domain
%     eD2{i} = zeros(L(end-i),1);      
%     eDr2{i} = zeros(len,1);          
%     delays(i) = 2^i-1;              % Level delay for synthesis
%     w{i} = zeros(L(end-i),1);       % Subband adaptive filter coefficient, initialize to zeros
%     w2{i} = zeros(L(end-i),1);
% end 
% w{i} = zeros(L(end-i),2);           % Last level has 2 columns, cD and cA
% eD{i} = zeros(1,2);                 % Last level has 2 columns, cD and cA
% U.tmp = zeros(len,1);
% Y.tmp = zeros(len,1);
% U.Z = zeros(2,1);
% Y.Z = zeros(2,1);
% 
% u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  
% y = zeros(len,1);                 % Tapped-delay line of desired response (Analysis FB)
% 
% ITER = length(un);
% en = zeros(1,ITER);               % Initialize error sequence to zero
% 
% 
% 
% % % ONLY FOR TESTING PURPOSE
% % t=0:0.001:1;
% % dn=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  
% 
% % Testing freezed filters
% % w{1} = zeros(L(end-1),2);
% % w{1}(1,:) = 1;
% 
% for n = 1:ITER    
%     u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
%     y = [dn(n); y(1:end-1)];        % Desired response vector        
% 
%     % Analysis Bank
%     U.tmp = u;
%     Y.tmp = y;
%     U.tmp2 = u;
%     Y.tmp2 = y;
%     for i = 1:level
%         if mod(n,2^i) == 0
%             U.Z = H'*U.tmp;
%             U.cD{i} = [U.Z(2); U.cD{i}(1:end-1)]; 
%             U.cA{i} = [U.Z(1); U.cA{i}(1:end-1)];
%             U.tmp = U.cA{i}(1:len);
%             U.Z2 = H'*U.tmp2;
%             U.cD2{i} = [U.Z2(2); U.cD2{i}(1:end-1)]; 
%             U.cA2{i} = [U.Z2(1); U.cA2{i}(1:end-1)];
%             U.tmp2 = U.cD2{i}(1:len);
%             
%             Y.Z = H'*Y.tmp;
%             Y.cD{i} = [Y.Z(2); Y.cD{i}(1:end-1)]; 
%             Y.cA{i} = [Y.Z(1); Y.cA{i}(1:end-1)];
%             Y.tmp = Y.cA{i}(1:len);
%             Y.Z2 = H'*Y.tmp2;
%             Y.cD2{i} = [Y.Z2(2); Y.cD2{i}(1:end-1)]; 
%             Y.cA2{i} = [Y.Z2(1); Y.cA2{i}(1:end-1)];
%             Y.tmp2 = Y.cD2{i}(1:len);        
%             
%             if i == level
%                 eD{i} = Y.Z' - sum(([U.cA{i}, U.cD{i}]).*w{i});
%                 eD2{i} = Y.Z2' - sum(([U.cA2{i}, U.cD2{i}]).*w2{i});                
% 
%                 if n >= AdaptStart(i)
%                     w{i} = w{i} + [U.cA{i},U.cD{i}].*(eD{i}./(sum([U.cA{i},U.cD{i}].*[U.cA{i},U.cD{i}])+alpha))*mu; 
%                     w2{i} = w2{i} + [U.cA2{i},U.cD2{i}].*(eD2{i}./(sum([U.cA2{i},U.cD2{i}].*[U.cA2{i},U.cD2{i}])+alpha))*mu; 
%                 end 
%             end           
%             S.iter{i} = S.iter{i} + 1;                
%         end
%     end    
% 
% 
%     % Synthesis Bank
%     for i = level:-1:1
%         if i == level
%             if mod(n,2^i) == 0
%                 eDr{i} = F*eD{i}' + eDr{i};
%                 eDr2{i} = F*eD2{i}' + eDr2{i};
%             end
%         else
%             if mod(n,2^i) == 0                
%                 eDr{i} = F*[eDr{i+1}(1); eDr2{i+1}(1)] + eDr{i};
%                 eDr{i+1} = [eDr{i+1}(2:end); 0];
%                 eDr2{i+1} = [eDr2{i+1}(2:end); 0];
%             end            
%         end
%     end   
%     en(n) = eDr{i}(1);
%     eDr{i} = [eDr{i}(2:end); 0];           
% end
% 
% en = en(1:ITER);
% S.coeffs = [w; w2];
% 
% 
% end


