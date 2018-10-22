function S = Volterra_Init(M, mu, level, wtype)

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
S.kernel_length  = M;              % Kernel dimensions
S.step          = mu;             % Step saize 
S.iter          = cell(1,level);  % Itertion count per level
S.levels        = level(1);          % DWT levels 
S.wtype         = wtype;          % Filter type
S.alpha         = 1e-6;           % Small positive constant

disp("only first level parameter is taken for entire structure for now, must be generalized")

[low_d,high_d,low_r,high_r] = wfilters(wtype);
S.analysis = [low_d', high_d'];     % Analysis filters
S.synthesis = [low_r', high_r'];    % Synthesis filters

M1 = M(1);
M2 = M(2)*(M(2)+1)/2;
S.length = [M1, M2];

S.AdaptStart = max(S.kernel_length);


end
