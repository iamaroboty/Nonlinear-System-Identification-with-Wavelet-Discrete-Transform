function S = QMFInit(M, mu, level, wtype)

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
S.step          = mu;             % Step size 
S.iter          = cell(1,level);  % Iteration count per level
S.levels        = level;          % DWT levels 
S.wtype         = wtype;          % Filter type
S.alpha         = 1e-6;           % Small positive constant

% [low_d,high_d,low_r,high_r] = wfilters(wtype);
% S.analysis = [qmf(low_d'), qmf(high_d')];     % Analysis filters
% S.synthesis = [qmf(low_r'), qmf(high_r')];    % Synthesis filters

% Create filter from QMF banks
db = dbwavf(wtype);
qmfdb = qmf(db); 
S.analysis = [qmfdb; db]'*sqrt(2); 
S.synthesis = flip(S.analysis);

% Subbands filter lenghts
lf = length(S.analysis);
L = [M; zeros(level,1)];
for i= 1:level
    L = [floor((L(1)+lf-1)/2); L(1:end-1)];
end
S.L = [L(1); L]';

for i=1:level
    S.AdaptStart(i) = 2^i*L(end-i);
end

% Oldshit
% L = zeros(level,1);
% for i= 1:level
%     L = [M/(2^i); L(1:end-1)];
%     S.iter{i} = 0;
% end
% S.L = [L(1); L; M]';            % Wavelet decomposition lenghts [cAn cDn cDn-1 ... cD1]



end
