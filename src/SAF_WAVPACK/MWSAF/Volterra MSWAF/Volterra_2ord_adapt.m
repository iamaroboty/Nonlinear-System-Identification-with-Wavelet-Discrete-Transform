
function [en,S] = Volterra_2ord_adapt(un,dn,S)
% Wavelet-Decomposition Subband Adaptive Filter (WAF)                 
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal

M = S.length;                     % kernel lengths  (array) 
mu = S.step;                      % Step Size here is an array 

AdaptStart = S.AdaptStart;        % Transient
alpha = S.alpha;                  % Small constant (1e-6)

% filter bank parameters
H = S.analysis;                   % Analysis filter bank
level = S.levels(1);                 % Wavelet Levels

display("only first level parameter is taken for entire structure for now, must be generalized")

H = create_multilevel_bank(H, level);
F = flip(H);

S.analysis_bank = H;
S.synthesis_bank = F;

[len, ~] = size(H);               % Wavelet filter length
L = S.L;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

U1 = zeros(S.length(1),2^level); % output of analysis filters for 1st order signal 
U2 = zeros(S.length(2),2^level); % output of analysis filters for 2nd order signal 

a1 = zeros(len,1);
a2 = zeros(len+S.length(2),1);
d = zeros(len,1); 
A1 = zeros(len,2^level);  
A2 = zeros(len,2^level);

z = zeros(len,1);

w1 = zeros(S.length(1),1);
w2 = zeros(S.lenght(2),1);

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero

%  t=0:0.001:1;
%  un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  
%  dn = un;
%  tot_delay = (2^level - 1)*(len-1) +1 ;

	
for n = 1:ITER
    
    d = [dn(n); d(1:end-1)];                       % Update tapped-delay line of d(n)
    
    a1 = [un(n); a(1:end-1)];                       % Update tapped-delay line of u(n)
    A1 = [a, A(:,1:end-1)];                         % Update buffer
    
    a2 = [un(n); a(1:end-1)];                       % Update tapped-delay line of u(n)
    A2 = [a, A(:,1:end-1)];                         % Update buffer
   
       
    if (mod(n,2^level)==0)                               % Tap-weight adaptation at decimated rate
        U1 = (H'*A1)';                              % Partitioning u(n) 
        U1_2 = U1(1:end-2^level,:);
        U1_tot = [U1', U1_2']';                           % Subband data matrix
        
        U2 = (H'*A2)';                              % Partitioning u(n) 
        U2_2 = U2(1:end-2^level,:);
        U2_tot = [U2', U2_2']';                           % Subband data matrix
        
        dD = H'*d;                                 % Partitioning d(n) 
        
        eD = dD - U1_tot'*w1- U2_tot'*w2;                         % Error estimation
        
        if n >= AdaptStart
            w1 = w1 + U1_tot*(eD./(sum(U1_tot.*U1_tot)+alpha)')*mu(1); % Tap-weight adaptation
            w2 = w2 + U2_tot*(eD./(sum(U2_tot.*U2_tot)+alpha)')*mu(2);
            
        end
        z = F*eD + z;                                       
        en(n-2^level+1:n) = z(1:2^level); 
        z = [z(2^level+1:end); zeros(2^level,1)]; 
    end
                          
end

S.coeffs = w;
end


    


