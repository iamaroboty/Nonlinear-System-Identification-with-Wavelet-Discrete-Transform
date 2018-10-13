function S = MWSAFinit(M, mu, level, wtype, Q)

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
S.iter          = 0;              % Itertion count per level
S.levels        = level;          % DWT levels 
S.wtype         = wtype;          % Filter type
S.alpha         = 1e-6;           % Small positive constant

if Q == 0
    [low_d,high_d,low_r,high_r] = wfilters(wtype);
    H = [low_d', high_d'];     % Analysis filters
    F = [low_r', high_r'];    % Synthesis filters
    S.analysis = H;
    S.synthesis = F;
%     E = reshape(low_d,2,length(low_d)/2);         % Each row represent a polyphase component
elseif Q == 1
    % Create filter from QMF banks
    db = dbwavf(wtype);
    qmfdb = qmf(db); 
    H = [qmfdb; db]'*sqrt(2); 
    S.analysis = H;
    S.synthesis = flip(H);%     E = reshape(db,2,length(db)/2);         % Each row represent a polyphase component
else
    error("Set Q either 1 or 0");
end

if level == 2
    Hi = upsample(H,2);
    Hi = [conv(Hi(:,1),H(:,1)), conv(Hi(:,2),H(:,1)), conv(Hi(:,1),H(:,2)), conv(Hi(:,2),H(:,2))];  
    if mod(length(Hi),2) ~= 0
        Hi = Hi(1:end-1,:);
    end
    S.analysis = Hi;
    S.synthesis = flip(Hi);
end



% Subbands filter lenghts
lf = length(S.analysis);
L = [M; zeros(level,1)];
for i= 1:level
    L = [floor((L(1)+lf-1)/2); L(1:end-1)];
end
S.L = [L(1); L]'./2;

% % Old
% L = zeros(level,1);
% for i= 1:level
%     L = [M/(2^i); L(1:end-1)];
%     S.iter{i} = 0;
% end
% S.L = [L(1); L; M]';            % Wavelet decomposition lenghts [cAn cDn cDn-1 ... cD1]

% for i=1:level
%     S.AdaptStart(i) = 2^i*S.L(end-i);
% end

S.AdaptStart = M+length(H);

% S.Polyphase = E;                % Polyphase representation of analysis filter bank

end

