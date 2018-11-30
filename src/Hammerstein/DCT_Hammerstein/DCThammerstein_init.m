function S = DCThammerstein_init(M, mu, order, alpha, N, L, D )

% DCThammerstein_init  Initialize Parameter Structure for DCT 
%                      Adaptive Hammerstein Filter
%
% Arguments: 
% M              Unknown sequence lenght (equaivalent adaptive filter
%                lenght)
% mu             Step size for LMS algorithm 
% order          taylor approx max order
% N              number of subbands
% L              lenght of analysis filters

if nargin < 7
    D = N ; % critical decfac
end

% Defalt filter bank: Pseudo-QMF CMFB

[hopt, ~] = opt_filter(L-1,N); % Generate a prototype lowpass filter
[H,F] = make_bank(hopt,N);           % Generate filter banks using cosine modulation
H = sqrt(N)*H';                      % Analysis section
F = sqrt(N)*F';                      % Synthesis section


% Assign structure fields
S.step          = mu;               % Step size
S.decfac        = D;                % Decimation factor
S.analysis      = H;                % Analysis filter bank
S.synthesis     = F;                % Synthesis filter bank
S.iter          = 0;                % Iteration count
S.alpha         = alpha;             % regularization positive constant
S.filt_length  = M;                  % lin_filt dimensions
S.AdaptStart    = L + max(M);       % Running effect of the analysis and adaptive 
                                    %   filter, minimum L + M
S.N_subbands = N ;         % n sub
S.order         = order;            % max power series hammersdtein order 

end

