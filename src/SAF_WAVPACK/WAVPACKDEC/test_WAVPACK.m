%% WAVELET PACKET DECOMPOSITION TEST FOR PERFECT RECONSTRUTION

addpath 'Common';             % Functions in Common folder
% clear all; close all

% Testing Signal

d = 256;        %Total signal length
t=0:0.001:10;
un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
un = un(1:d)';


%% wavpack parameters

mu = 0.3;                      % ignored here 
M = 256;                        % Length of unknown system response also ignored here
level = 2;                     % Levels of Wavelet decomposition
filters = 'db2';               % Set wavelet type
Ovr = 1;                       % Oversampling factor

% S = QMFInit(M,mu,level,filters);
S = SWAFinit(M, mu, level, filters); 

M = S.length;                     % Unknown system length (Equivalent adpative filter lenght)

F = S.analysis;                   % Analysis filter bank
H = S.synthesis;                  % Synthesis filter bank

% analysis and synthesis are used in reverse to obtain in U.Z a column
% vector with cD in the first position


[len, ~] = size(H);               % Wavelet filter length
level = S.levels;                 % Wavelet Levels
L = S.L.*Ovr;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

% Init Arrays
for i= 1:level     
       U.c{i} = zeros(L(end-i),2^(i));            
       eD{i} = zeros(L(end-i),2^(i-1));      % Error signa, transformed domain
       eDr{i} = zeros(len,2^(i-1));          % Error signal, time domain
       delays(i) = 2^i-1;                    % Level delay for synthesis
    
          
end 
w = zeros(L(end-i),2^i);           % Last level has 2 columns, cD and cA

w(1,1:end) = 1;                    % set filters to kronecker delta

eD{i} = zeros(1,2^i);              % Last level has 2 columns, cD and cA

pwr = w;
beta = 1./L(2:end-1);

u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


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
            
            U.Z = HH'*U.tmp; % column [cD ; cA] 
            
         
            [rows, cols] = size(U.Z);
            
            indx = 1;
            
            for col=1:cols
                for row=1:rows 
                    
                    U.c{i}(:,indx) = cat(1,U.Z(row,col), U.c{i}(1:end-1, indx)); %CA||CD
        
                indx=indx+1;
                end  
            end
            
            U.tmp = U.c{i}(1:len,:);    
             
            
            if i == level
                
                eD{i} = sum((U.c{i}).*w);        
        
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
               
                for col = 1:2:2^i-1 
                 
                    
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


