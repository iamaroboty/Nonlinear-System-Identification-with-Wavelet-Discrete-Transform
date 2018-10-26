  function [en,S] = SAFTERRA_adapt(un,dn,S)

% SAFadapt          Subband Adaptive Filter
%
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in MSAFinit.m
% en                History of error signal
%
% by Lee, Gan, and Kuo, 2008
% Subband Adaptive Filtering: Theory and Implementation
% Publisher: John Wiley and Sons, Ltd

[L,N] = size(S.analysis);
D = S.decfac;
mu = S.step;
alpha = S.alpha;
AdaptStart = S.AdaptStart;
H = S.analysis;
F = S.synthesis;

W = S.coeffs; 

% Adaptive subfilters
w1 = zeros(size(W{1})); 

U = zeros(size(w1));           % Tapped-delay lines of adaptive subfilters

for i = 1:size(S.coeffs{2},1)
    u2{i} = zeros(size(W{2},1),1);   % Nonlinear input delay line
    
    U2{i} = zeros(size(W{2},1), size(W{2},2));     
    w2{i} = zeros(size(W{2},1)-i+1, size(W{2},2));
    
end


H = H(:,1:N/2+1);             % Analysis filter (only the first N/2+1 are taken)
x = zeros(size(W{2},1),1);               % Tapped-delay line of input signal (Analysis FB)
y = zeros(L,1);               % Tapped-delay line of desired response (Analysis FB)
z = zeros(L,1);               % Tapped-delay line of error signal (Synthesis FB)

ITER = length(un);
en = zeros(1,ITER);

for n = 1:ITER
    x = [un(n); x(1:end-1)];  % Fullband input vector for band partitioning
    y = [dn(n); y(1:end-1)];  % Fullband desired response vector for band partitioning
    
    for i = 1:size(S.coeffs{2},1)
    u2{i} = [un(n)*x(i); u2{i}(1:end-1)]; %input vectors for second order term 
    end
    
    if (mod(n-1,D)==0)              % Tap-weight adaptation at lowest sampling rate
        
        U = [x(1:L)'*H; U(1:end-1,:)];   % Each colums hold a subband regression vector
        
        ker2 = 0; 
        norm = 0; 
        for i = 1:size(S.coeffs{2},1)
        U2{i} = [u2{i}(1:L)'*H; U2{i}(1:end-1,:)];  
        ker2 = ker2 + sum(U2{i}(1:size(w2{i},1),:).*w2{i});
        norm = norm + sum(U2{i}(1:size(w2{i},1),:).*conj(U2{i}(1:size(w2{i},1),:)));
        end
        
        dD = y'*H;                  % Row vector, each element is decimated desired signal
        eD = dD - sum(U.*w1) - ker2;
        if n >= AdaptStart
            w1 = w1 + conj(U)*diag(eD./(sum(U.*conj(U))+alpha))*mu(1);
            
            for i= 1:size(S.coeffs{2},1)
                w2{i} = w2{i} + conj(U2{i}(1:size(w2{i},1),:))*diag(eD./(norm +alpha))*mu(2);
             
            end
            
            
            S.iter = S.iter + 1;
        end
        eD = [eD, conj(fliplr(eD(2:end-1)))].'; 
        z = F*eD + z;                                       
    end
    en(n) = real(z(1)); 
    z = [z(2:end); 0];
                   
end

en = en(1:ITER);
S.coeffs{1} = w1;
S.coeffs{2} = w2; 
