function [en, S] = modwtmraHammerstein_adapt(un, dn, S)

w = S.coeffs;
M = S.length;
Iter = S.iter;
AdaptStart = S.AdaptStart;
palpha = S.palpha;
alpha = S.alpha;
wv = S.wavelet;

d = zeros(M,1);
u = zeros(M,1);
ITER = length(un);  
yn = zeros(1,ITER); 
en = zeros(1,ITER); 

J = floor(log2(M)) +1;
p = ones(J,1);
w = zeros(M,1);


for n = 1:ITER
    u = [un(n); u(1:end-1)];
    d = [dn(n); d(1:end-1)];
    
    if (mod(n,M) == 0)
        wtecg = modwt(u, wv);
        U = modwtmra(wtecg, wv);
        
        est = sum(U.*W);
        
        e = d - est';
        
        if n >= AdaptStart
            W = W + U*diag(eD./(sum(U.*U,2)+alpha))*mu;
        end
    end

end

