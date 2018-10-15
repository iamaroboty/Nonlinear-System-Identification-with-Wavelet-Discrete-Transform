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



% petraglia structure extended analysis filters
% 4 band structure 
H_1lvl = cat(2, conv(H(:,1), upsample(H(:,1),2)), conv(H(:,1), upsample(H(:,2),2)), conv(H(:,2), upsample(H(:,1),2)), conv(H(:,2), upsample(H(:,2),2)));
% H0 , H1, H2, H3
H_af = cat(2, conv(H_1lvl(:,1), H_1lvl(:,1)), conv(H_1lvl(:,1), H_1lvl(:,2)), conv(H_1lvl(:,2), H_1lvl(:,2)),  conv(H_1lvl(:,3), H_1lvl(:,2)),  conv(H_1lvl(:,3), H_1lvl(:,3)),conv(H_1lvl(:,3), H_1lvl(:,4)),  conv(H_1lvl(:,4), H_1lvl(:,4)) );
% H0H0 , H0H1, H1H1, H1H2, H2H2; 
F_1lvl = cat(2, conv(F(:,1), upsample(F(:,1),2)), conv(F(:,1), upsample(F(:,2),2)), conv(F(:,2), upsample(F(:,1),2)), conv(F(:,2), upsample(F(:,2),2)));

% high pass , band pass, low pass 


[len_af, ~] = size(H_af);               % Wavelet filter length
[len, ~] = size(H); 


%% petraglia aliasing free structure adaptation



% analysis and synthesis are used in reverse to obtain in U.Z a column
% vector with cD in the first position


% Init Arrays
for i= 1:level
     
       U.c{i} = zeros(L(end-i),7);    
    
       Y.c{i} = zeros(L(end-i),2^(i));   % nb is different !!!    
    
       eD{i} = zeros(L(end-i),2^(i-1));      % Error signa, transformed domain
       eDr{i} = zeros(len,2^(i-1));          % Error signal, time domain
       delays(i) = 2^i-1;                    % Level delay for synthesis
    
      
end 
w = zeros(L(end-i),4);           % Last level has 2 columns, cD and cA
eD{i} = zeros(1,2^i/2);              % Last level has 2 columns, cD and cA

pwr = w;
beta = 1./L(2:end-1);

u = zeros(len_af,1);                 % Tapped-delay line of input signal (Analysis FB)  
y = zeros(len,1);                 % Tapped-delay line of desired response (Analysis FB)

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


% % ONLY FOR TESTING PURPOSE
%  t=0:0.001:1;
%  un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  
%  ITER = length(un);
%  w(1,1:end) = 1;

ytap = zeros(len_af-len,1); 

for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    y = [ytap(end); y(1:end-1)];        % Desired response vector
    
    ytap = [dn(n); ytap(1:end-1)]; 
    
    % Analysis Bank
    U.tmp = u;
    Y.tmp= y;
    
    for i = 1:level
        if mod(n,2^i) == 0
    
            U.Z = H_af'*U.tmp;
            Y.Z = H'*Y.tmp;  % normal
            
            [rows, cols] = size(U.Z);
            
            indx = 1;
            
            for col=1:cols
                for row=1:rows 
                    
                    U.c{i}(:,indx) = cat(1,U.Z(row,col), U.c{i}(1:end-1, indx));
                    
              
                    indx=indx+1;
                end  
            end
            
            indx = 1;
            [rows, cols] = size(Y.Z);
            
             for col=1:cols
                for row=1:rows 
                    
                    Y.c{i}(:,indx) = cat(1,Y.Z(row,col), Y.c{i}(1:end-1, indx));
                    
              
                    indx=indx+1;
                end  
            end
            
                   
            
            %U.tmp = U.c{i}(1:len,:);    
            %Y.tmp = Y.c{i}(1:len,:);  
            
            if i == level
                
                
                directcD = sum(U.c{i}(:,1).*w(:,1)); 
                directcA = sum(U.c{i}(:,3).*w(:,2)); 
                
                
                crosscD = sum(U.c{i}(:,2).*w(:,2)); 
                crosscA = sum(U.c{i}(:,2).*w(:,1));
                
                
                filtered =  [directcD+crosscD; directcA+crosscA]'; 

      
                eD{i} = Y.Z' - filtered; 
                
%                  for col=1:cols
%                    for row=1:rows
%                       eD{i}(indx) = Y.Z(row,col) - filtered(indx);
%                        
%                
%                       indx=indx+1;
%                    end  
%                  end



                                
                if n >= AdaptStart(i)
%                     pwr{i} = beta(i)*pwr{i}+ (1-beta(i))*([U.cA{i},U.cD{i}].*[U.cA{i},U.cD{i}]);
%                     w{i} = w{i} + mu*[U.cA{i},U.cD{i}].*eD{i}./((sum(pwr{i})+alpha)); 
                    
                   % w = w + mu*U.c{i}.*eD{i}./(sum(U.c{i}.*U.c{i})+alpha); 
                    
                    w(:,1) = w(:,1) + mu*(U.c{i}(:,1).*eD{i}(1)./(sum(U.c{i}(:,1).*U.c{i}(:,1))+alpha) +  U.c{i}(:,2).*eD{i}(2)./(sum(U.c{i}(:,2).*U.c{i}(:,2))+alpha) ); 
                    w(:,2) = w(:,2) + mu*(U.c{i}(:,3).*eD{i}(2)./(sum(U.c{i}(:,3).*U.c{i}(:,3))+alpha) +  U.c{i}(:,2).*eD{i}(1)./(sum(U.c{i}(:,2).*U.c{i}(:,2))+alpha) ); 
%                     w(:,1) = 0.707.*w(:,1) ;
%                     w(:,2) = 0.707.*w(:,2) ;
%                                                                             
                                                
                                                
                end                     
            end           
            S.iter{i} = S.iter{i} + 1;                
        end
    end    

      % Synthesis Bank
    for i = level:-1:1
            
        if i == level
            if mod(n,2^i) == 0
                indx = 1; 
               
                for col = 1:2:size(eD{i},2)-1 
                 
                    
                eDr{i}(:,indx) = F*eD{i}(1,col:col+1)' + eDr{i}(:,indx);
                indx = indx +1;
                
                end
                
                
            end
        else
            if mod(n,2^i) == 0     
                
                indx = 1; 
                
                for col = 1:2:2^i-1
                eDr{i}(:,indx) = F*eDr{i+1}(1,col:col+1)' + eDr{i}(:,indx);
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
S.coeffs = w;
end


