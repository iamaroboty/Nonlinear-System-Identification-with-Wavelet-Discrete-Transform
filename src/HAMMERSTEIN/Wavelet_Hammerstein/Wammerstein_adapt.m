function [en,S] = Wammerstein_adapt(un,dn,S)
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
H = S.analysis;                   % Analysis filter bank
level = S.levels;                 % Wavelet Levels

H = create_multilevel_bank(H, level);
F = flip(H);

S.analysis_bank = H;
S.synthesis_bank = F;

[len, ~] = size(H);               % Wavelet filter length
L = S.L;                     % Wavelet decomposition Length, sufilter length [cAn cDn cDn-1 ... cD1 M]

order = S.order; 
xp = zeros(order, 1);  
w = zeros(M,1);
p = zeros(order, 1); % non linearity coeffs vector 
p(1) = 1; 
X = zeros(M, order); 

C = zeros(M,2^level,order);    
% C1 = zeros(2^level,2^level,order);                  
a = zeros(len,1); d = zeros(len,1); A = zeros(len,2^level);    % Analysis filtering
z = zeros(len,1);
tot_delay = (2^level - 1) * (len - 1);

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero

	
for n = 1:ITER
    
    d = [dn(n); d(1:end-1)];                       % Update tapped-delay line of d(n)
    a = [un(n); a(1:end-1)];                       % Update tapped-delay line of u(n)
    A = [a, A(:,1:end-1)];                         % Update buffer    
    
    if (mod(n,2^level)==0)                         % Tap-weight adaptation at decimated rate
        
        for i = 1:order
            C1 = (H'*A.^i)';
            C2 = C(1:end-2^level,:,i);
%             C = cat(1,C1,C2);
            C(:,:,i) = [C1', C2']';               % Taylor series expansion
            D(:,:,i) = C(:,:,i).*p(i);            % FIR filtering with p
        end
        
        U = sum(D,3);                             % FIR filtering with p
        
        dD = H'*d;                                 % Partitioning d(n) 
        eD = dD - U'*w;                            % Error estimation
        
        if n >= AdaptStart
            Ctmp = permute(C, [2 1 3]);         % generalization of traspose [C] = [M x 2^level x order]
            norm = 0;
            for i = 1:order
                tmp(i,:) = eD'*Ctmp(:,:,i);
                norm = norm + Ctmp(:,:,i).*Ctmp(:,:,i); %NOT SURE
            end
            
            p = p + (mu(1)*tmp*w)/((norm*w)'*(norm*w) + alpha); %NOT SURE
            w = w + U*(eD./(sum(U.*U)+alpha)')*mu(2); 
            
            S.iter = S.iter + 1;
        end
        
        z = F*eD + z;                                       
%         en(n-2^level+1:n) = z(1:2^level); 
%         z = [z(2^level+1:end); zeros(2^level,1)]; 
    end
    en(n) = z(1);
    z = [z(2:end); 0];                     
end

S.coeffs{1} = w;
S.coeffs{2} = p; 
end


    


