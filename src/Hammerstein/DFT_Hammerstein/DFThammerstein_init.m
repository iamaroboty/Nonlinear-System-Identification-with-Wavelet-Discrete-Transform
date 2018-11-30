function S = DFThammerstein_init(M, mu, order, alpha, N,D,L)

% DFThammerstein_init  Initialize Parameter Structure for DFT 
%                      Adaptive Hammerstein Filter
%
% Arguments: 
% M              Unknown sequence lenght (equaivalent adaptive filter
%                lenght)
% mu             Step size for LMS algorithm 
% order          taylor approx max order
% N              number of subbands
% D              decimation factor
% L              lenght of analysis filters


hopt = fir1(L-1,1/N);               % Generate prototype lowpass filter
[H,F] = make_bank_DFT(hopt,N);      % Generate filter banks using complex modulation
H = sqrt(D)*H;                      % Analysis section
F = sqrt(D)*F;                      % Synthesis section


% Assign structure fields
S.order = order;                    % taylor max order 
S.n_subbands = N;
S.step          = mu;               % Step size
S.decfac        = D;                % Decimation factor
S.analysis      = H;                % Analysis filter bank
S.synthesis     = F;                % Synthesis filter bank
S.iter          = 0;                % Iteration count
S.alpha         = alpha;             % Small positive constant 
S.coeffs{1}    = zeros(M/D,N/2+1);  % linear system filter len
S.coeffs{2} = zeros(order, 1);         % power series (taylor) len 
S.AdaptStart    = L + max(M);     % Running effect of the analysis and adaptive 
                              



end

