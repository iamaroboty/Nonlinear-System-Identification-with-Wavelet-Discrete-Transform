function S = SAFTERRA_Init(M,mu,N,D,L)

% SAFTERRA_Init       Initialize Parameter Structure for the Subband Adaptive Volterra Filter
%                    DFT Filter Bank Is Used by Default
%
% Arguments:
% M             Length of corresponding fullband filter
% mu            Step size
% N             Number of subbands
% D             Decimation factor, D=N (Critical decimation), D<N (oversampling)
%
% Subband Adaptive Filtering: Theory and Implementation
% by Lee, Gan, and Kuo, 2008
% Publisher: John Wiley and Sons, Ltd

% Defualt is DFT filter bank (Section 2.6)

hopt = fir1(L-1,1/N);               % Generate prototype lowpass filter
[H,F] = make_bank_DFT(hopt,N);      % Generate filter banks using complex modulation
H = sqrt(D)*H;                      % Analysis section
F = sqrt(D)*F;                      % Synthesis section

% Assign structure fields
S.step          = mu;               % Step size
S.decfac        = D;                % Decimation factor
S.analysis      = H;                % Analysis filter bank
S.synthesis     = F;                % Synthesis filter bank
S.iter          = 0;                % Iteration count
S.alpha         = 1e-6;             % Small positive constant
S.kernel_length  = M;              % Kernel dimensions
S.coeffs{1}     = zeros(M(1)/D,N/2+1);
S.coeffs{2}   =  zeros(max(M(2)/D,L),N/2+1); % NOT SURE 
S.AdaptStart    = L + max(M);            % Running effect of the analysis and adaptive 
                                    %   filter, minimum L + M

end
