%% WAVELET PACKET DECOMPOSITION TEST FOR PERFECT RECONSTRUTION

addpath 'Common';             % Functions in Common folder
clear all; close all

% Testing Signal

d = 256;        %Total signal length
t=0:0.001:10;
un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
un = sin(2*pi*t*10);
% un = ones(d,1);
% un = zeros(d,1); un(1) = 1;
un = un(1:d);
Ovr = 1; 


%% wavpack parameters

mu = 0.1;                      % ignored here 
M = 256;                        % Length of unknown system response also ignored here
level = 2;
filters = 'db1';               % Set wavelet type


% S = QMFInit(M,mu,level,filters);
S = SWAFinit(M, mu, level, filters); 

% db = dbwavf(filters);
% qmfdb = qmf(db); 
% H = [qmfdb; db]'; 
% F = flip(H);

F = S.analysis;                   % Analysis filter bank
H = S.synthesis;                  % Synthesis filter bank 
%% petraglia aliasing free structure adaptation

% filters for the aliasing free bank 

%upsampled filters


% %check for two layers
% check_H_extd = cat(2, conv(H(:,1), upsample(H(:,1),2)), conv(H(:,1), upsample(H(:,2),2)), conv(H(:,2), upsample(H(:,1),2)), conv(H(:,2), upsample(H(:,2),2)) ); 
% % H0, H1, H2, H3, H4
% % Hi = upsample(H,2);
% % Hi = [conv(Hi(:,1),H(:,1)), conv(Hi(:,2),H(:,1)), conv(Hi(:,1),H(:,2)), conv(Hi(:,2),H(:,2))];  
% % if mod(length(Hi),2) ~= 0
% %     Hi = Hi(1:end-1,:);
% % end
% % S.analysis = Hi;
% % S.synthesis = flip(Hi);
% % 
% 
% check_H_extd = check_H_extd(1:end-1,:);
% 
% 
% check_H_af = cat(2, conv(check_H_extd(:,1), check_H_extd(:,1)), conv(check_H_extd(:,1), check_H_extd(:,2)), ...
%                 conv(check_H_extd(:,2), check_H_extd(:,2)),  conv(check_H_extd(:,2), check_H_extd(:,3)), ...
%                   conv(check_H_extd(:,3), check_H_extd(:,3)),  conv(check_H_extd(:,3), check_H_extd(:,4)), ...
%                    conv(check_H_extd(:,4), check_H_extd(:,4))   ); 



Hi = zeros(2^(level-1)*size(H,1), 2^(level));

indx = 1;
for i = 1:size(H,2)
for j=1:level-1
    
   up{i,j} = upsample(H(:,i), 2^(j)); 
   
   
end
end

% up = upsample(H,2);
%outer product
H_tmp = H; 
for i=1:size(up,2)
H_tmp = outer_conv(H_tmp, up(:,i));
end

Hi = H_tmp; 

% if mod(length(Hi),2) ~= 0
%     Hi = Hi(1:end-1,:);
% end

HH = create_multilevel_bank(H,level);

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

%H_af = cat(2, conv(H(:,1), H(:,1)), conv(H(:,1), H(:,2)), conv(H(:,2), H(:,2))); 

F = flip(Hi);
eq = [6.125, 4.8125, 6.5625, 4.375];
% F = F./eq;


% figure;
% for i = 1:7
% plot(10*log10(abs(fft(H_af(:,i),512))), 'LineWidth',2); hold on;
% end
% legend('H0H0', 'H0H1', 'H1H1' , 'H1H2', 'H2H2', 'H2H3', 'H3H3');
% title('Petraglia Structure');
% axis([-inf 256 -40 inf])
% 
% figure;
% for i = 1:4
% plot(10*log10(abs(fft(Hi(:,i),512))), 'LineWidth',2); hold on;
% end
% legend('H0', 'H1', 'H2','H3');
% title('2 level filterbank (db1)');
% axis([-inf 256 -20 inf])

% analysis and synthesis are used in reverse to obtain in U.Z a column
% vector with cD in the first position


H_af = H_af(1:end-2,:);

[len_af, ~] = size(H_af);               % Wavelet filter length
[len, ~] = size(Hi); 

L = S.L.*Ovr;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

% Init Arrays

% everything is brought to the first level
    
U_c = zeros(L(end-level),2^(level+1)-1);            
eDr = zeros(len,1);          % Error signal, time domain
tmp = zeros(len_af,1); 
delay = 1;                    % Level delay for synthesis
z = zeros(len,1);
           
w = zeros(L(end-level),2^level);           % Last level has 2 columns, cD and cA
w(1,:) = 1/sqrt(4) ;

eD = zeros(1,2^level);              % Last level has 2 columns, cD and cA

% pwr = w;
% beta = 1./L(2:end-1);

u = zeros(len_af,1);                 % Tapped-delay line of input signal (Analysis FB)  

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero

for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'

    % Analysis Bank
    U.tmp = u;
   
        if (mod(n,2^level) == 0)            
            
            U.Z = H_af'*U.tmp; % column [cD ; cA] 
%             U.Z = Hi'*U.tmp;
         
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
            crossH3H2G2 = sum(U_c(:,6).*w(:,2));
            
            summed = [directH0H0+crossH1H0G1; directH1H1+crossH1H0G0+crossH2H1G2;...
                        directH2H2+crossH2H1G1+crossH3H2G3; ...
                        directH3H3+crossH3H2G2];
            

            eD =  (summed) ;    
%             eD = U.Z;
            
           z = F*eD + z;                                       
           en(n-2^level+1:n) = z(1:2^level); 
           z = [z(2^level+1:end); zeros(2^level,1)]; 
            
        end
%         en(n) = eDr(1);
%         eDr = [eDr(2:end); 0];  
        
        
      
end
en = en(1:ITER);


%% check for perfect reconstruction

tot_delay = len_af-1;

stem(en(tot_delay:end));
hold on;
stem(un); 

