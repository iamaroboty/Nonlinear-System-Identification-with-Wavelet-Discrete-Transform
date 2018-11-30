function [en,S] = DFThammerstein_adapt(un,dn,S)
% DFT Subband Adaptive Hammerstein Filter (WAF)                 
% 
% Arguments:
% un                Input signal
% dn                Desired signal
% S                 Adptive filter parameters as defined in WSAFinit.m
% en                History of error signal


mu = S.step;                      % Step Size
AdaptStart = S.AdaptStart;        % Transient
N_subbands = S.n_subbands;        % number of subbands
decfac = S.decfac;                     % decimation factor
alpha = S.alpha;                  % Small constant (1e-6)
H = S.analysis;                   % Analysis filter bank
F = S.synthesis;                  % Synthesis filter bank


[len, ~] = size(H);          % DFT analysis filter length

H = H(:,1:N_subbands/2+1); % Analysis filter (only the first N/2+1 are taken)


order = S.order; 

w = S.coeffs{1};
p = S.coeffs{2}; % non linearity coeffs vector 
p(1) = 1; 

C = zeros(size(w,1),size(w,2), order);                 
a = zeros(len,1); d = zeros(len,1); A = zeros(len, size(w,2));  
z = zeros(len,1);

ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero


% [L,N] = size(S.analysis);
% D = S.decfac;
% mu = S.step;
% alpha = S.alpha;
% AdaptStart = S.AdaptStart;
% H = S.analysis;
% F = S.synthesis;
% W = S.coeffs{1};                 % Adaptive subfilters
% U = zeros(size(w));           % Tapped-delay lines of adaptive subfilters
% H = H(:,1:N/2+1);             % Analysis filter (only the first N/2+1 are taken)
% x = zeros(L,1);               % Tapped-delay line of input signal (Analysis FB)
% y = zeros(L,1);               % Tapped-delay line of desired response (Analysis FB)
% z = zeros(L,1);               % Tapped-delay line of error signal (Synthesis FB)
% 
% ITER = length(un);
% en = zeros(1,ITER);
% 
% for n = 1:ITER
%     x = [un(n); x(1:end-1)];  % Fullband input vector for band partitioning
%     y = [dn(n); y(1:end-1)];  % Fullband desired response vector for band partitioning
%     
%     if (mod(n-1,D)==0)              % Tap-weight adaptation at lowest sampling rate
%         U = [x'*H; U(1:end-1,:)];   % Each colums hold a subband regression vector
%         dD = y'*H;                  % Row vector, each element is decimated desired signal
%         eD = dD - sum(U.*W);
%         if n >= AdaptStart
%             W = W + conj(U)*diag(eD./(sum(U.*conj(U))+alpha))*mu(2);
%             S.iter = S.iter + 1;
%         end
%         eD = [eD, conj(fliplr(eD(2:end-1)))].'; 
%         z = F*eD + z;                                       
%     end
%     en(n) = real(z(1)); 
%     z = [z(2:end); 0];
%                    
% end

% en = en(1:ITER);
% S.coeffs = W;
% t=0:0.001:1;
% dn=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  
%  dn = un;
	


for n = 1:ITER
    
    d = [dn(n); d(1:end-1)];                       % Update tapped-delay line of d(n)
    a = [un(n); a(1:end-1)];                       % Update tapped-delay line of u(n)
    A = [a, A(:,1:end-1)];                         % Update buffer    
    
    if (mod(n-1,decfac)==0)                         % Tap-weight adaptation at decimated rate
        
        for i = 1:order
            C1 = (H'*A.^i)';    %(H'*A.^i)';            % C1 contains all taps for different hammerstein-taylor orders
            C2 = C(1:end-size(w,2),:,i);      
            C(:,:,i) = [C1', C2']';     % Taylor series expansion
            D(:,:,i) = C(:,:,i).*p(i);  % FIR filtering with p polynomial 
        end
        
        U = sum(D,3);                             % FIR filtering with p

        
        dD = d'*H;                                 % Partitioning d(n) 
        eD = dD - sum(U.*w);                            % Error estimation
        
        if n >= AdaptStart
            Ctmp = permute(C, [2 1 3]);   % generalization of traspose [C] = [M x N_subbands x order]
            %norm = 0;
            for i = 1:order
                %tmp =  sum((eD*Ctmp(:,:,i)*w).*conj(eD*Ctmp(:,:,i)*w)); %eD*sum(Ctmp(:,:,i)*w)';
                %norm = sum((Ctmp(:,:,i)*w).*(conj(Ctmp(:,:,i)*w))); 
                %p(i) = p(i) + sum(conj(sum(Ctmp(:,:,i)*w))*diag(eD./((sum(Ctmp(:,:,i)*w)*sum(Ctmp(:,:,i)*w)') + alpha)))*mu(1);%(mu(1)*tmp)/( sum(norm.*norm) + alpha);
                p(i) = p(i) + sum(conj(sum(Ctmp(:,:,i)*w))*diag(eD./(sum((sum(Ctmp(:,:,i)*w)).*conj(sum(Ctmp(:,:,i)*w)))+alpha))*mu(1)); 
                %norm = norm + sum(sum((Ctmp(:,:,i)*w).*(conj(Ctmp(:,:,i)*w)))); 
            end
            %(mu(1)*en(n)*X'*w)/((X'*w)'*(X'*w) + alpha)
             
            w = w +  conj(U)*diag(eD./(sum(U.*conj(U))+alpha))*mu(2); 
         
            S.iter = S.iter + 1;
        end
        
        eD = [eD, conj(fliplr(eD(2:end-1)))].';
        z = F*eD + z;                                       

    end
    en(n) = real(z(1));
    z = [z(2:end); 0];                     
end

S.coeffs{1} = w;
S.coeffs{2} = p; 
end


    


