function S = MSAFTERRA_Init(M,mu,N,L)

if nargin > 4
    if (size(H,2)~=N)|(size(F,2)~=N)
        error('Columns of H (%d) or F (%d) not match with N = %d',size(H,2),size(F,2),N);
    end
else

% Defualt filter bank: Pseudo-QMF CMFB

    [hopt,passedge] = opt_filter(L-1,N); % Generate a prototype lowpass filter
    [H,F] = make_bank(hopt,N);           % Generate filter banks using cosine modulation
    H = sqrt(N)*H';                      % Analysis section
    F = sqrt(N)*F';                      % Synthesis section
end

% Assign structure fields
S.step          = mu;               % Step size
S.decfac        = N;                % Decimation factor
S.analysis      = H;                % Analysis filter bank
S.synthesis     = F;                % Synthesis filter bank
S.iter          = 0;                % Iteration count
S.alpha         = 1e-6;             % Small positive constant
S.kernel_length  = M;              % Kernel dimensions
S.AdaptStart    = L + max(M);            % Running effect of the analysis and adaptive 
                                    %   filter, minimum L + M

end
