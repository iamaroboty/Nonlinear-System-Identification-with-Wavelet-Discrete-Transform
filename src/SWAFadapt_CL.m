function [en,S] = SWAFadapt_CL(un,dn,S)
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
H = S.analysis;                   % Analysis filter bank
F = S.synthesis;                  % Synthesis filter bank
[len, ~] = size(H);               % Wavelet filter length
level = S.levels;                 % Wavelet Levels
L = S.L;                          % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

% Init Arrays
for i= 1:level
    U.cD{i} = zeros(L(end-i),1);    
    U.cA{i} = zeros(L(end-i),1); 
    
    Ut.cD{i} = zeros(L(end-i),1);    
    Ut.cA{i} = zeros(L(end-i),1);
    
    E.cD{i} = zeros(L(end-i),1);
    E.cA{i} = zeros(L(end-i),1);
    
    eD{i} = zeros(L(end-i),1);  % Error signa, transformed domain
    YD{i} = zeros(L(end-i),1); 
    
    eDr{i} = zeros(len,1);          % Error signal, time domain
    YDr{i} = zeros(len,1);
    delays(i) = 2^i-1;              % Level delay for synthesis
    w{i} = zeros(L(end-i),1);       % Subband adaptive filter coefficient, initialize to zeros 
    %w{i}(1) = 1;
end 
w{i} = zeros(L(end-i),2);           % Last level has 2 columns, cD and cA
%w{i}(1,1) = 1;
%w{i}(1,2) = 1;
eD{i} = zeros(1,2);                 % Last level has 2 columns, cD and cA
YD{i} = zeros(1,2); 
U.tmp = zeros(len,1);
Ut.tmp = zeros(len,1);
Y.tmp = zeros(len,1);
U.Z = zeros(2,1);
Ut.Z = zeros(2,1);
Y.Z = zeros(2,1);

u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  
e = zeros(len,1); 
yn = zeros(len,1); % Tapped-delay line of desired response (Analysis FB)

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


% % ONLY FOR TESTING PURPOSE
% t=0:0.001:1;
% un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  

% Testing freezed filters
% w{1} = zeros(L(end-1),2);
% w{1}(1,:) = 1;


tot_delay = (2^level - 1)*(len-1) +1 ;

yn_tap = zeros(tot_delay,1); 
un_tap = zeros(tot_delay+len,1); 

dn_tap = zeros(tot_delay,1); 

U_cDtap = cell(level,tot_delay); 
U_cAtap = cell(level,tot_delay); 

%dn = horzcat([dn(1), dn(1), dn(1), dn(1), dn(1)] , dn); 
%un = horzcat([0.1, 0.1, 0.1, 0.1, 0.1] ,  dn); 

en_tap = zeros(len,1); 

for n = 1:ITER  
    
    % current value goes on HEAD 
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    en_tap = [dn_tap(tot_delay)-yn_tap(1); en_tap(1:end-1)];  %y = [dn(n); y(1:end-1)];        % time error vector  
    
    dn_tap = [dn(n);  dn_tap(1:end-1)];   
    un_tap = [un(n);  un_tap(1:end-1)];   
    
    % Analysis Bank
    U.tmp = u;
    Ut.tmp = un_tap(end-len+1:end);
    E.tmp = en_tap;
    for i = 1:level
        if mod(n,2^i) == 0
            U.Z = H'*U.tmp;
            U.cD{i} = [U.Z(2); U.cD{i}(1:end-1)]; 
            U.cA{i} = [U.Z(1); U.cA{i}(1:end-1)];
            U.tmp = U.cA{i}(1:len);
            
            Ut.Z = H'*Ut.tmp;
            Ut.cD{i} = [Ut.Z(2); Ut.cD{i}(1:end-1)]; 
            Ut.cA{i} = [Ut.Z(1); Ut.cA{i}(1:end-1)];
            Ut.tmp = Ut.cA{i}(1:len);
            
            
            E.Z = H'*E.tmp;
            E.cD{i} = [E.Z(2); E.cD{i}(1:end-1)]; 
            E.cA{i} = [E.Z(1); E.cA{i}(1:end-1)];
            E.tmp = E.cA{i}(1:len);
            
%             
%             for index= tot_delay:-1:2
%                 U_cDtap{i,index} = U_cDtap{i,index-1}; % 2*4 
%                 U_cAtap{i,index} = U_cAtap{i,index-1}; %2*8 
%             end
%          
%             U_cDtap{i,1} = U.cD{i};
%             U_cAtap{i,1} = U.cA{i};
%             
            
            
            if i == level
                %eD{i} =  [E.cA{i}(end), E.cD{i}(end)];
                YD{i} = sum(([U.cA{i}, U.cD{i}]).*w{i}); % dot product 

                if n >= AdaptStart+ tot_delay
                    w{i} = w{i} + [Ut.cA{i},Ut.cD{i}].*(E.Z'./(sum([Ut.cA{i},Ut.cD{i}].*[Ut.cA{i},Ut.cD{i}])+alpha))*mu; 
                end 
            else
                YD{i} = [YD{i}(2:end); U.cD{i}'*w{i}];  
   

                if n >= AdaptStart + tot_delay
                    w{i} = w{i} + (mu*E.Z(2)/(Ut.cD{i}'*Ut.cD{i} + alpha))*Ut.cD{i};
                end
            end           
            S.iter{i} = S.iter{i} + 1;                
        end
    end    


    % Synthesis Bank
    for i = level:-1:1
        if i == level
            if mod(n,2^i) == 0
                YDr{i} = F*YD{i}' + YDr{i};
            end
        else
            if mod(n,2^i) == 0                
                YDr{i} = F*[YDr{i+1}(1); YD{i}(end-(len-1)*delays(end-i))] + YDr{i};
                YDr{i+1} = [YDr{i+1}(2:end); 0];
            end            
        end
    end   
    yn_tap = [YDr{i}(1); yn_tap(1:end-1)];
    YDr{i} = [YDr{i}(2:end); 0]; % upsample
    en(n) = en_tap(1); 
    
    
    
end


en = en(1:ITER); 
S.coeffs = w;
end

% %% Full packet 2 layer : TESTING
% % Init Arrays
% for i= 1:level
%     U.cD{i} = zeros(L(end-i),1);    
%     U.cA{i} = zeros(L(end-i),1);    
%     Y.cD{i} = zeros(L(end-i),1);
%     Y.cA{i} = zeros(L(end-i),1);
%     eD{i} = zeros(L(end-i),1);      % Error signa, transformed domain
%     eDr{i} = zeros(len,1);          % Error signal, time domain
%     delays(i) = 2^i-1;              % Level delay for synthesis
%     w{i} = zeros(L(end-i),1);       % Subband adaptive filter coefficient, initialize to zeros
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
% % un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  
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
%                 if n >= AdaptStart
%                     w{i} = w{i} + [U.cA{i},U.cD{i}].*(eD{i}./(sum([U.cA{i},U.cD{i}].*[U.cA{i},U.cD{i}])+alpha))*mu; 
% <<<<<<< HEAD
%                     w2{i} = w2{i} + [U.cA2{i},U.cD2{i}].*(eD2{i}./(sum([U.cA2{i},U.cD2{i}].*[U.cA2{i},U.cD2{i}])+alpha))*mu; 
% =======
%                     
% >>>>>>> a504b831d16ba8d3704244303f3af62d8ac996d8
%                 end 
%             else
%                 eD{i} = [eD{i}(2:end); Y.cD{i}(1) - U.cD{i}'*w{i}]; 
% 
%                 if n >= AdaptStart
%                     w{i} = w{i} + (mu*eD{i}(end)/(U.cD{i}'*U.cD{i} + alpha))*U.cD{i};
%                 end
%             end
%            
%             S.iter{i} = S.iter{i} + 1;
%                 
%         end
%     end    
% 
% 
%     % Synthesis Bank
%     for i = level:-1:1
%         if i == level
%             if mod(n,2^i) == 0
%                 eDr{i} = F*eD{i}' + eDr{i};
%             end
%         else
%             if mod(n,2^i) == 0                
%                 eDr{i} = F*[eDr{i+1}(1); eD{i}(end-(len-1)*delays(end-i))] + eDr{i};
%                 eDr{i+1} = [eDr{i+1}(2:end); 0];
%             end            
%         end
%     end   
%     en(n) = eDr{i}(1);
%     eDr{i} = [eDr{i}(2:end); 0];           
% end
% 
% en = en(1:ITER);
% S.coeffs = w;
% 
% 
% end


