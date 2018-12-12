function S = modwtmraHammerstein_init(M, mu, wv, palpha, alpha)



% Assign structure fields only after error checking is complete

S.coeffs        = zeros(M,1);     % Coefficients of FIR filter 
S.length        = M;            % Length of adaptive filter
S.step          = mu;           % Step size of unconstrained FDAF 
S.iter          = 0;          % Iteration count
S.AdaptStart    = M;            % Running effect of adaptive filter
S.palpha        = palpha;        % Constant to update the power in each frequency bin
S.alpha = alpha; 
% S.level = lev;
S.wavelet = wv;
                              
end

