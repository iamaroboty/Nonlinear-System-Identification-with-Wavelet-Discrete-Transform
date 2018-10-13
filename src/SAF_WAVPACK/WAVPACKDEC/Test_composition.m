%% TEST PACKET DECOMPOSITION
%% WAVELET PACKET DECOMPOSITION TEST FOR PERFECT RECONSTRUTION

addpath 'Common';             % Functions in Common folder
clear all; close all

% Testing Signal

d = 256;        %Total signal length
t=0:0.001:10;
un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
un = un(1:d)';


mu = 0.3;                      % ignored here 
M = 256;                        % Length of unknown system response also ignored here
level = 2;                     % Levels of Wavelet decomposition
filters = 'db1';               % Set wavelet type
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


%% NORMAL 
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

%% Noble identity
%% init H
Hi = upsample(H,2);
Hi = [conv(Hi(:,1),H(:,1)), conv(Hi(:,2),H(:,1)), conv(Hi(:,1),H(:,2)), conv(Hi(:,2),H(:,2))];    
Hi = Hi(1:end-1,:);
H = Hi;
F = flip(Hi);
%%

[len, ~] = size(H);               % Wavelet filter length
level = 1;                 % Wavelet Levels
L = S.L.*Ovr./2;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

% Init Arrays
for i= 1:level     
       UU.c{i} = zeros(L(end-i),4^(i));        
       Y.c{i} = zeros(L(end-i),4^(i));        
       eD{i} = zeros(L(end-i),4^(i-1));      % Error signa, transformed domain
       eDr{i} = zeros(len,4^(i-1));          % Error signal, time domain
       delays(i) = 4^i-1;                    % Level delay for synthesis    
       %w{i} = zeros(L(end-i),1);       % Subband adaptive filter coefficient, initialize to zeros    
end 
w = zeros(L(end-i),4^i);   
w(1,1:end) = 1;
w_cross = zeros(L(end-i),6);
eD{i} = zeros(1,4^i);     


eDcross = zeros(1,6);
cross_filtered = zeros(1,4^i);

% pwr = w;
% beta = 1./L(2:end-1);

u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  
y = zeros(len,1);                 % Tapped-delay line of desired response (Analysis FB)

ITER = length(un);
enn = zeros(1,ITER);               % Initialize error sequence to zero


% % ONLY FOR TESTING PURPOSE
%  t=0:0.001:1;
%  un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  
%  ITER = length(un);
% Testing freezed filters
% w{1} = zeros(L(end-1),2);
% 
%  w(1,1:end) = 1;

for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'

    % Analysis Bank
    UU.tmp = u;
    
    for i = 1:level
        if mod(n,4^i/(Ovr)) == 0
            if (i==1 && Ovr == 2)
                HH = H./sqrt(2);
            else
                HH = H;
            end
            
            UU.Z = HH'*UU.tmp;
            
            [rows, cols] = size(UU.Z);
            
            indx = 1;
            
            for col=1:cols
                for row=1:rows                     
                    UU.c{i}(:,indx) = cat(1,UU.Z(row,col), UU.c{i}(1:end-1, indx));
              
                    indx=indx+1;
                end  
            end
            
            UU.tmp = UU.c{i}(1:len,:);    
            Y.tmp = Y.c{i}(1:len,:);  
            
            if i == level                
                filtered = sum((UU.c{i}).*w);     
%                 
%                 crossUcD1 = sum(UU.c{i}(:,2).*w_cross(:,1));  % UU.cD1
%                 crossUcA1 = sum(UU.c{i}(:,1).*w_cross(:,2)) +  sum(UU.c{i}(:,3).*w_cross(:,3));  
%                 crossUcD2 = sum(UU.c{i}(:,2).*w_cross(:,4)) +  sum(UU.c{i}(:,4).*w_cross(:,5)); 
%                 crossUcA2 = sum(UU.c{i}(:,3).*w_cross(:,6));
%                 
%                 cross_filtered = [ crossUcD1; crossUcA1;  crossUcD2; crossUcA2];
%                 
%                 cross_Uc = cat(2, UU.c{i}(:,2), UU.c{i}(:,1), UU.c{i}(:,3), UU.c{i}(:,2), UU.c{i}(:,4), UU.c{i}(:,3));
                
                indx=1;
                
                 for col=1:cols
                   for row=1:rows
%                       eD{i}(indx) = Y.Z(row,col) - (filtered(indx)+cross_filtered(indx));    
                      eD{i}(indx) = (filtered(indx));
                                      
                      indx=indx+1;
                   end  
                 end
                 
                 eDcross = [eD{i}(1), eD{i}(2), eD{i}(2), eD{i}(3), eD{i}(3), eD{i}(4)];
                                
%                 if n >= 256
% %                     pwr{i} = beta(i)*pwr{i}+ (1-beta(i))*([UU.cA{i},UU.cD{i}].*[UU.cA{i},UU.cD{i}]);
% %                     w{i} = w{i} + mu*[UU.cA{i},UU.cD{i}].*eD{i}./((sum(pwr{i})+alpha)); 
%                     
%                     w = w + mu*UU.c{i}.*eD{i}./(sum(UU.c{i}.*UU.c{i})+alpha); 
%                     w_cross = w_cross + mu*cross_Uc.*eDcross./(sum(cross_Uc.*cross_Uc)+alpha);     
%                                                 
%                 end                     
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
            if mod(n,4^i/(Ovr)) == 0
                    eDr{i} = FF*eD{i}' + eDr{i};                       
            end
        else
            if mod(n,4^i/(Ovr)) == 0     
                
                indx = 1;                 
                for col = 1:4:4^i-1
                    eDr{i}(:,indx) = FF*eDr{i+1}(1,col:col+1)' + eDr{i}(:,indx);
                    indx = indx +1;                
                end
                
                eDr{i+1} = [eDr{i+1}(2:end,:); zeros(1,size(eDr{i+1},2))];
            end            
        end
    end   
    enn(n) = eDr{i}(1);
    eDr{i} = [eDr{i}(2:end); 0];           
end

%% check for perfect reconstruction

tot_delay = (2^level - 1)*(len-1) +1 ;

stem(enn(tot_delay:end));
hold on;
stem(un);





