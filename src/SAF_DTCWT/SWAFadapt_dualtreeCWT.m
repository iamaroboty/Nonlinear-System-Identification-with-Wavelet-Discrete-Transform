function [en,S] = SWAFadapt_dualtreeCWT(un,dn,S)
% SWAFadapt         Wavelet-transformed Subband Adaptive Filter (WAF)                 
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
%H = S.analysis;                   % Analysis filter bank
%F = S.synthesis;                  % Synthesis filter bank

level = S.levels;                 % Wavelet Levels
L = S.L;                          % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]




df = dtfilters('dtf1');

fdf = df{1}; 
df = df{2}; 

[len, ~] = size(fdf{1});
%[len2, ~] = size(df{1}{1});



% Init Arrays
for i= 1:level
    
    if level == 1
    H1{i} =  fdf{1};
    H2{i} = fdf{2};
    F1{1} = flipud(fdf{1});
    F2{1} = flipud(fdf{2});
    
    
    else
    
    H1{i} = df{1};
    H2{i} = df{2};
    F1{i} = flipud(df{1});
    F2{i} = flipud(df{2});
    end
    
    U1.cD{i} = zeros(L(end-i),1);    
    U1.cA{i} = zeros(L(end-i),1);   
    U2.cD{i} = zeros(L(end-i),1);    
    U2.cA{i} = zeros(L(end-i),1);    
    Y1.cD{i} = zeros(L(end-i),1);
    Y1.cA{i} = zeros(L(end-i),1);
    Y2.cD{i} = zeros(L(end-i),1);
    Y2.cA{i} = zeros(L(end-i),1);
    
    eD1{i} = zeros(L(end-i),1);      % Error signa, transformed domain
    
   
    eDr1{i} = zeros(len,1);          % Error signal, time domain
   
    
    
    eD2{i} = zeros(L(end-i),1);      % Error signa, transformed domain
    eDr2{i} = zeros(len,1);          % Error signal, time domain
    
    delays(i) = 2^i-1;              % Level delay for synthesis
    w1{i} = zeros(L(end-i),1);       % Subband adaptive filter coefficient, initialize to zeros    
    
    w2{i} = zeros(L(end-i),1);
    
end 
w1{i} = zeros(L(end-i),2);           % Last level has 2 columns, cD and cA
eD1{i} = zeros(1,2);                 % Last level has 2 columns, cD and cA
U1.tmp = zeros(len,1);
Y1.tmp = zeros(len,1);
U1.Z = zeros(2,1);
Y1.Z = zeros(2,1);

w2{i} = zeros(L(end-i),2);           % Last level has 2 columns, cD and cA
eD2{i} = zeros(1,2);                 % Last level has 2 columns, cD and cA
U2.tmp = zeros(len,1);
Y2.tmp = zeros(len,1);
U2.Z = zeros(2,1);
Y2.Z = zeros(2,1);

%pwr = w;
%beta = 1./L(2:end-1);

u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  
y = zeros(len,1);                 % Tapped-delay line of desired response (Analysis FB)

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


% % ONLY FOR TESTING PURPOSE
% t=0:0.001:1;
% un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  

% Testing freezed filters
% w{1} = zeros(L(end-1),2);
% w{1}(1,:) = 1;

for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    y = [dn(n); y(1:end-1)];        % Desired response vector        

    % Analysis Bank
    U1.tmp = u;
    Y1.tmp = y;
    
    U2.tmp = u;
    Y2.tmp = y;
    
    for i = 1:level
        if mod(n,2^i) == 0
    
            U1.Z = H1{i}'*U1.tmp;
            U1.cD{i} = [U1.Z(2); U1.cD{i}(1:end-1)]; 
            U1.cA{i} = [U1.Z(1); U1.cA{i}(1:end-1)];
            U1.tmp = U1.cA{i}(1:len);
            
            U2.Z = H2{i}'*U2.tmp;
            U2.cD{i} = [U2.Z(2); U2.cD{i}(1:end-1)]; 
            U2.cA{i} = [U2.Z(1); U2.cA{i}(1:end-1)];
            U2.tmp = U2.cA{i}(1:len);
            
            
            Y1.Z = H1{i}'*Y1.tmp;
            Y1.cD{i} = [Y1.Z(2); Y1.cD{i}(1:end-1)]; 
            Y1.cA{i} = [Y1.Z(1); Y1.cA{i}(1:end-1)];
            Y1.tmp = Y1.cA{i}(1:len);
            
            Y2.Z = H2{i}'*Y2.tmp;
            Y2.cD{i} = [Y2.Z(2); Y2.cD{i}(1:end-1)]; 
            Y2.cA{i} = [Y2.Z(1); Y2.cA{i}(1:end-1)];
            Y2.tmp = Y2.cA{i}(1:len);
            
            
            if i == level
                eD1{i} = Y1.Z' -  sum(([U1.cA{i}, U1.cD{i}]).*w1{i}) ;
                
                %Y = complex(Y1.Z', Y2.Z');
                
                %UcA = complex(U1.cA{i}, U2.cA{i}) ;
                %UcD = complex(U1.cD{i} ,U2.cD{i}) ;
                
                %eD =  Y - sum(([UcA, UcD]).*w1{i}); 
                
                %eD1{i} = real(eD); 
                %eD2{i} = imag(eD); %Y2.Z' 
                
                eD2{i} = Y2.Z'- sum(([U2.cA{i}, U2.cD{i}]).*w2{i});
                
                if n >= AdaptStart(i)
%                     pwr{i} = beta(i)*pwr{i}+ (1-beta(i))*([U.cA{i},U.cD{i}].*[U.cA{i},U.cD{i}]);
%                    w{i} = w{i} + mu*[U.cA{i},U.cD{i}].*eD{i}./((sum(pwr{i})+alpha)); 
                w1{i} = w1{i} + mu*[U1.cA{i},U1.cD{i}].*eD1{i}./(sum([U1.cA{i},U1.cD{i}].*[U1.cA{i},U1.cD{i}])+alpha) ; 
               %w1{i} = w1{i} + mu*[UcA,UcD].*eD./(sum([UcA,UcD].*[UcA,UcD])+alpha) ; 
                w2{i} = w2{i} + mu*[U2.cA{i},U2.cD{i}].*eD2{i}./(sum([U2.cA{i},U2.cD{i}].*[U2.cA{i},U2.cD{i}])+alpha); 
                    
                end 
            else
                eD1{i} = [eD1{i}(2:end); Y1.cD{i}(1) - U1.cD{i}'*w1{i}]; 
                
                
                eD2{i} = [eD2{i}(2:end); Y2.cD{i}(1) - U2.cD{i}'*w2{i}]; 
                
                
                if n >= AdaptStart(i)
%                     pwr{i} = beta(i)*pwr{i}+ (1-beta(i))*(U.cD{i}.*U.cD{i});
%                     w{i} = w{i} + (mu*eD{i}(end)/((sum(pwr{i}) + alpha)))*U.cD{i};
                    w1{i} = w1{i} + (mu*eD1{i}(end)/(U1.cD{i}'*U1.cD{i} + alpha))*U1.cD{i};
                    w2{i} = w2{i} + (mu*eD2{i}(end)/(U2.cD{i}'*U2.cD{i} + alpha))*U2.cD{i};
                    
                end
            end           
            S.iter{i} = S.iter{i} + 1;                
        end
    end    

    % Synthesis Bank
    for i = level:-1:1     
        if i == level
            if mod(n,2^i) == 0
                eDr1{i} = F1{i}*eD1{i}' + eDr1{i};
                eDr2{i} = F2{i}*eD2{i}' + eDr2{i};
                
            end
        else
            if mod(n,2^i) == 0                
                eDr1{i} = F1{i}*[eDr1{i+1}(1); eD1{i}(end-(len-1)*delays(end-i))] + eDr1{i};
                eDr1{i+1} = [eDr1{i+1}(2:end); 0];
                
                eDr2{i} = F2{i}*[eDr2{i+1}(1); eD2{i}(end-(len-1)*delays(end-i))] + eDr2{i};
                eDr2{i+1} = [eDr2{i+1}(2:end); 0];
            end            
        end
    end   
    en(n) = 0.5*(eDr1{i}(1) + eDr2{i}(1));
    eDr1{i} = [eDr1{i}(2:end); 0];       
    eDr2{i} = [eDr2{i}(2:end); 0];  
end

en = en(1:ITER);
S.coeffs = [w1, w2];
end

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


