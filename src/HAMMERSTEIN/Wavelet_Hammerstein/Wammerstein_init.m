function S = Wammerstein_init(M, mu, level, wtype, order, alpha)

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
S.length        = M;              % Unknown system length (Equivalent adpative filter length)
% S.AdaptStart    = M;              % Transient
S.step          = mu;             % Step saize 
S.iter          = 0;  % Itertion count per level
S.levels        = level;          % DWT levels 
S.wtype         = wtype;          % Filter type
S.alpha         = alpha;           % Small positive constant
S.AdaptStart    = M;
S.order         = order;


[low_d,high_d,low_r,high_r] = wfilters(wtype);
S.analysis = [low_d', high_d'];     % Analysis filters
S.synthesis = [low_r', high_r'];    % Synthesis filters

% Subbands filter lenghts
lf = length(low_d);
L = [M; zeros(level,1)];
for i= 1:level
    L = [floor((L(1)+lf-1)/2); L(1:end-1)];
end
S.L = [L(1); L]';

% Oldshit
% L = zeros(level,1);
% for i= 1:level
%     L = [M/(2^i); L(1:end-1)];
%     S.iter{i} = 0;
% end
% S.L = [L(1); L; M]';            % Wavelet decomposition lenghts [cAn cDn cDn-1 ... cD1]

end

