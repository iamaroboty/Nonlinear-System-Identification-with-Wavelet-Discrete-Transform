function S = Hammerstein_NLMS_init(M, mu, leaks)

% SWAFinit          Initialize Parameter Structure for Wavelet Transform
%                   Adaptive Filter
%
% Arguments: 
% M              Unknown sequence lenght (equaivalent adaptive filter
%                lenght)
% mu             Step size for LMS algorithm 
% level          Levels of the Wavelet decomposition (dwt)
% wtype          Mother Wavelet filter type

% Assign structure fields
S.filters_lengths  = M;              % Kernel dimensions
S.step          = mu;             % Step saize 
S.alpha         = 1e-6;           % Small positive constant
S.leaks = leaks; 
S.iter = 0; 
S.AdaptStart = max(M);

end
