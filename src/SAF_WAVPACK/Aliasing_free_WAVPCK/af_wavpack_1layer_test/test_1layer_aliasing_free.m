%% WAVELET PACKET DECOMPOSITION TEST FOR PERFECT RECONSTRUTION

addpath 'Common';             % Functions in Common folder
clear all; close all

% Testing Signal

d = 256;        %Total signal length
t=0:0.001:10;
un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
% un = ones(1,d);
un = un(1:d);
Ovr = 1; 


%% wavpack parameters

mu = 0.1;                      % ignored here 
M = 256;                        % Length of unknown system response also ignored here
level = 1;                     % Levels of Wavelet decomposition
filters = 'db8';               % Set wavelet type


%S = QMFInit(M,mu,level,filters);
S = SWAFinit(M, mu, level, filters); 

M = S.length;                     % Unknown system length (Equivalent adpative filter lenght)

F = S.analysis;                   % Analysis filter bank
H = S.synthesis;                  % Synthesis filter bank 


%% petraglia aliasing free structure adaptation

% filters for the aliasing free bank 
H_af = cat(2, conv(H(:,1), H(:,1)), conv(H(:,1), H(:,2)), conv(H(:,2), H(:,2))); 
%F_af = cat(2, conv(F(:,1), F(:,1)), conv(F(:,2), F(:,2)), conv(F(:,2), F(:,2))); 

% analysis and synthesis are used in reverse to obtain in U.Z a column
% vector with cD in the first position


[len_af, ~] = size(H_af);               % Wavelet filter length
[len, ~] = size(H); 

level = S.levels;                 % Wavelet Levels
L = S.L.*Ovr;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

% Init Arrays
for i= 1:level     
       U.c{i} = zeros(L(end-i),2^(i)+1);            
       eD{i} = zeros(L(end-i),2^(i-1));      % Error signa, transformed domain
       eDr{i} = zeros(len,2^(i-1));          % Error signal, time domain
       delays(i) = 2^i-1;                    % Level delay for synthesis              
end 
w = zeros(L(end-i),2^i);           % Last level has 2 columns, cD and cA

w(1,:) = 1/sqrt(2);                   % set filters to kronecker delta


eD{i} = zeros(1,2^i/2);              % Last level has 2 columns, cD and cA

pwr = w;
beta = 1./L(2:end-1);

u = zeros(len_af,1);                 % Tapped-delay line of input signal (Analysis FB)  

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero

tot_delay = (2^level - 1)*(len-1) +1 ;

%%
for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'

    % Analysis Bank
    U.tmp = u;
    
    for i = 1:level
        if mod(n,2^i/(Ovr)) == 0
            if (i==1 && Ovr == 2)
                HH = H./sqrt(2);
            else
                HH = H;
            end
            
            U.Z = H_af'*U.tmp; % column [cD ; cA] 
            
         
            [rows, cols] = size(U.Z);
            
            indx = 1;
            
            for col=1:cols
                for row=1:rows 
                    
                    U.c{i}(:,indx) = cat(1,U.Z(row,col), U.c{i}(1:end-1, indx)); %CD||CA
        
                indx=indx+1;
                end  
            end
            
            U.tmp = U.c{i}(1:len,:);    
             
            
            if i == level
                
                directcD = sum(U.c{i}(:,1).*w(:,1)); 
                directcA = sum(U.c{i}(:,3).*w(:,2)); 
                
                
                crosscD = sum(U.c{i}(:,2).*w(:,2)); 
                crosscA = sum(U.c{i}(:,2).*w(:,1));
                
%                 direct = sum(U.c{i}(:,1:2:end).*w(:,1:2:end)); 
%                 
%                 cross = sum(U.c{i}(:,2:2:end).*w(:,2:2:end)); 
                
                eD{i} = [directcD+crosscD; directcA+crosscA]';  
               
                
        
            end   
            S.iter{i} = S.iter{i} + 1;                
        end
    end    


    % Synthesis Bank
    for i = level:-1:1
            if (i==1 && Ovr == 2)
                FF = F./sqrt(2);
            else
                FF = F;
            end
        if i == level
            if mod(n,2^i/(Ovr)) == 0
                indx = 1; 
               
                for col = 1:2:size(eD{i},2)-1 
                 
                    
                eDr{i}(:,indx) = FF*eD{i}(1,col:col+1)' + eDr{i}(:,indx);
                indx = indx +1;
                
                end
                
                
            end
        else
            if mod(n,2^i/(Ovr)) == 0     
                
                indx = 1; 
                
                for col = 1:2:2^i-1
                eDr{i}(:,indx) = FF*eDr{i+1}(1,col:col+1)' + eDr{i}(:,indx);
                indx = indx +1;
                
                end
                eDr{i+1} = [eDr{i+1}(2:end,:); zeros(1,size(eDr{i+1},2))];
            end            
        end
    end   
    en(n) = eDr{i}(1);
    eDr{i} = [eDr{i}(2:end); 0];           
end

en = en(1:ITER);


%% check for perfect reconstruction

tot_delay = (2^level - 1)*(len-1) +1 ;

stem(en(tot_delay:end));
hold on;
stem(un); 

% %%
% U.c = zeros(L(end-level),2^(level+1)-1); 
% z = zeros(len,1);
% 
% for n = 1:ITER    
%     u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
% 
%     % Analysis Bank
%     U.tmp = u;
% %     U.tmp = u(1:len);
%     
%         if (mod(n,2^level) == 0)            
%             
%             U.Z = H_af'*U.tmp; % column [cD ; cA] 
% %             U.Z = Hi'*U.tmp;
%          
%             [rows, cols] = size(U.Z);
%             
%             indx = 1;
%             
%             for col=1:cols
%                 for row=1:rows 
%                     
%                     U.c(:,indx) = cat(1,U.Z(row,col), U.c(1:end-1, indx)); %CD||CA
%         
%                 indx=indx+1;
%                 end  
%             end
%             
%                 directcD = sum(U.c(:,1).*w(:,1)); 
%                 directcA = sum(U.c(:,3).*w(:,2)); 
%                 
%                 
%                 crosscD = sum(U.c(:,2).*w(:,2)); 
%                 crosscA = sum(U.c(:,2).*w(:,1));
%                 
%                 
%                 eD = [directcD+crosscD; directcA+crosscA];    
%             
%            z = F*eD + z;                                       
%            en(n-2^level+1:n) = z(1:2^level); 
%            z = [z(2^level+1:end); zeros(2^level,1)]; 
%             
%         end
% %         en(n) = eDr(1);
% %         eDr = [eDr(2:end); 0];  
%         
%         
%       
% end
% en = en(1:ITER);
