function [en,S] = MWAFadapt(un,dn,S)

% WSAFadapt         Wavelet-transformed Subband Adaptive Filter (WAF)
%                   Transformation with an orthogonal matrix W.                    
%
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal

M = length(S.coeffs);
mu = S.step;                     % Step Size
beta = S.beta;                   % Forgettig factor
AdaptStart = S.AdaptStart;
W = S.W;                         % Transform Matrix
lev = S.levels;                  % Wavelet Levels
wtype = S.wtype;                 % Wavelet family   
w = S.coeffs;                    % Adaptive filtering
b = S.unknownsys;

L = zeros(lev,1);                   % Each subband lenght [cAn cDn cDn-1 ... cD1]

u = zeros(M,1);
y = zeros(M,1);
z = zeros(M,1);

ITER = length(un);
en = zeros(1,ITER);                 % Initialize error sequence to zero

% wcD1 = zeros(M/2, 1);
% wcD2 = zeros(M/4, 1);
% wcD3 = zeros(M/8, 1);
% wcD4 = zeros(M/16, 1);
% wcA = zeros(M/16, 1);
% ecD1 = zeros(M/2, 1);
% ecD2 = zeros(M/4, 1);
% ecD3 = zeros(M/8, 1);
% ecD4 = zeros(M/16, 1);
% ecA = zeros(M/16, 1);
% pwr_cD1 = zeros(M/2, 1);
% pwr_cD2 = zeros(M/4, 1);
% pwr_cA = zeros(M/4, 1);

% Calculate Filter length
[low_d,~,~,~] = wfilters(wtype);
lf = length(low_d);
L = [M; zeros(lev,1)];
for i= 1:lev
    L = [floor((L(1)+lf-1)/2); L(1:end-1)];
end
L = [L(1); L]';

wcD1 = zeros(L(5), 1);
wcD2 = zeros(L(4), 1);
wcD3 = zeros(L(3), 1);
wcD4 = zeros(L(2), 1);
wcA = zeros(L(1), 1);
ecD1 = zeros(L(5), 1);
ecD2 = zeros(L(4), 1);
ecD3 = zeros(L(3), 1);
ecD4 = zeros(L(2), 1);
ecA = zeros(L(1), 1);

% % %little help
% wcD1 = [1; zeros(M/2 -1,1)];
% wcD2 = [1; zeros(M/4 -1,1)];
% wcA = [1; zeros(M/4 -1,1)];
% 
t=0:0.001:1;
un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
dwtmode('sym')
[UU,L] = wavedec(flip(un(1:256)), lev, wtype);
% UU = W*flip(un(1:256))';
for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
    y = [dn(n); y(1:end-1)];        % Desired response vector        
       
    if mod(n,2) == 0  
%         U = W*u;                        % Transformed (DWT) input vector [Mx1]
%         Y = W*y;                        % Transformed (DWT) y vector [Mx1]
        [U, L] = wavedec(u, lev, wtype);
        [Y, L] = wavedec(y, lev, wtype);

        % Decomposition
        [UcD1, UcD2, UcD3, UcD4] = detcoef(U, L, [1, 2, 3, 4]);
        [YcD1, YcD2, YcD3, YcD4] = detcoef(Y, L, [1, 2, 3, 4]);
%         UcA = U(1:M/(2^lev));
%         YcA = Y(1:M/(2^lev));        
        UcA = appcoef(U, L, wtype);
        YcA = appcoef(Y, L, wtype);
                 
        ecD1 = [YcD1(1) - UcD1'*wcD1; ecD1(1:end-1)];
        if n >= AdaptStart
        %                 pwr_cD1 = (1-1/(M/2))* pwr_cD1+(1/(M/2))*(UcD1.*UcD1);
        %                 wcD1 = wcD1 + (mu*ecD1(1))./sqrt(pwr_cD1+0.00001).*UcD1;  
            wcD1 = wcD1 + (mu*ecD1(1)/(UcD1'*UcD1 + 0.00001))*UcD1;
        end

        if mod(n,4) == 0           
            ecD2 = [YcD2(1) - UcD2'*wcD2; ecD2(1:end-1)];
            if n >= AdaptStart
                %                 pwr_cD1 = (1-1/(M/2))* pwr_cD1+(1/(M/2))*(UcD1.*UcD1);
                %                 wcD1 = wcD1 + (mu*ecD1(1))./sqrt(pwr_cD1+0.00001).*UcD1;  
                wcD2 = wcD2 + (mu*ecD2(1)/(UcD2'*UcD2 + 0.00001))*UcD2;
            end

            if mod(n,8) == 0           
                ecD3 = [YcD3(1) - UcD3'*wcD3; ecD3(1:end-1)];
                if n >= AdaptStart
                    %                 pwr_cD1 = (1-1/(M/2))* pwr_cD1+(1/(M/2))*(UcD1.*UcD1);
                    %                 wcD1 = wcD1 + (mu*ecD1(1))./sqrt(pwr_cD1+0.00001).*UcD1;  
                    wcD3 = wcD3 + (mu*ecD3(1)/(UcD3'*UcD3 + 0.00001))*UcD3;
                end

                if mod(n,16) == 0
                    ecD4 = [YcD4(1) - UcD4'*wcD4; ecD4(1:end-1)];
                    ecA = [YcA(1) - UcA'*wcA; ecA(1:end-1)];
                    if n >= AdaptStart
                        %                     pwr_cD2 = (1-1/(M/4))* pwr_cD2+(1/(M/4))*(UcD2.*UcD2);
                        %                     pwr_cA = (1-1/(M/4))* pwr_cA+(1/(M/4))*(UcA.*UcA);        
                        %                     wcD2 = wcD2 + (mu*ecD2(1))./sqrt(pwr_cD2+0.00001).*UcD2;
                        %                     wcA = wcA + (mu*ecA(1))./sqrt(pwr_cA+0.00001).*UcA;
                        wcD4 = wcD4 + (mu*ecD4(1)/(UcD4'*UcD4 + 0.00001))*UcD4;
                        wcA = wcA + (mu*ecA(1)/(UcA'*UcA + 0.00001))*UcA;
                        S.iter = S.iter + 2;   
                    end
                    ew = [ecA; ecD4; ecD3; ecD2; ecD1];                                                         
%                     z = W'*ew; 
                    z = waverec(ew, L, wtype);
                    en(n-16+1:n) = flip(z(1:16));
                end  
            end 
        end
    end

%     en(n) = z(1);
    
    if mod(n,5000)== 0
        plot(10*log10(en(1:n).^2));
        xlabel('Number of iteration'); 
        ylabel('Live MSE error (dB)');linkdata on    %Live plotting      
    end
    
end


en = en(1:ITER);
S.coeffs = w;
end


    
