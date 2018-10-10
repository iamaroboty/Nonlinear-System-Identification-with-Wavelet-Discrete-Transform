function S = DSWAFinit(M, mu, level, wtype, J, Q)

% SWAFinit          Initialize Parameter Structure for Wavelet Transform
%                   Adaptive Filter
%
% Arguments: 
% M              Unknown sequence lenght (equaivalent adaptive filter
%                lenght)
% mu             Step size for LMS algorithm 
% level          Levels of the Wavelet decomposition (dwt)
% wtype          Mother Wavelet filter type
% Q              Flag for QMF type filter bank (e.g. filters obtainted with
%                qmf Matlab function)

% Assign structure fields
S.length        = M;              % Unknown system length (Equivalent adpative filter length)
S.step          = mu;             % Step saize 
S.iter          = cell(1,level);  % Itertion count per level
S.levels        = level;          % DWT levels 
S.wtype         = wtype;          % Filter type
S.alpha         = 1e-6;           % Small positive constant
S.UpdateRate    = round(M/J);     % J is typically in the range 1 to 8

if Q == 0
    [low_d,high_d,low_r,high_r] = wfilters(wtype);
    S.analysis = [low_d', high_d'];     % Analysis filters
    S.synthesis = [low_r', high_r'];    % Synthesis filters
%     E = reshape(low_d,2,length(low_d)/2);         % Each row represent a polyphase component
elseif Q == 1
    % Create filter from QMF banks
    db = dbwavf(wtype);
    qmfdb = qmf(db); 
    S.analysis = [qmfdb; db]'*sqrt(2); 
    S.synthesis = flip(S.analysis);
%     E = reshape(db,2,length(db)/2);         % Each row represent a polyphase component
else
    error("Set Q either 1 or 0");
end

% Subbands filter lenghts
% lf = length(S.analysis);
% L = [M; zeros(level,1)];
% for i= 1:level
%     L = [floor((L(1)+lf-1)/2); L(1:end-1)];
% end
% S.L = [L(1); L]';

% % Old
L = zeros(level,1);
for i= 1:level
    L = [M/(2^i); L(1:end-1)];
    S.iter{i} = 0;
end
S.L = [L(1); L; M]';            % Wavelet decomposition lenghts [cAn cDn cDn-1 ... cD1]

for i=1:level
    S.AdaptStart(i) = 2^i*S.L(end-i);
end

% S.Polyphase = E;                % Polyphase representation of analysis filter bank

end

