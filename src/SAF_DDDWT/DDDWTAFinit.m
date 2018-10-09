function S = DDDWTAFinit(M, mu, level, wtype)

% DDDWTAFinit      Initialize Parameter Structure for Double Density Wavelet Transform
%                   Adaptive Filter
%
% Arguments: 
% M              Unknown sequence lenght (equaivalent adaptive filter
%                lenght)
% mu             Step size for LMS algorithm 
% level          Levels of the Wavelet decomposition (dwt)
% filters          Mother Wavelet filter type

% Assign structure fields
S.length        = M;              % Unknown system length (Equivalent adpative filter length)
% S.AdaptStart    = M;              % Transient
S.step          = mu;             % Step saize 
S.iter          = cell(1,level);  % Itertion count per level
S.levels        = level;          % DWT levels 
S.filters         = wtype;          % Filter type
S.alpha         = 1e-6;           % Small positive constant

H = dtfilters(wtype);
F = flipud(H);  

S.analysis = H;
S.synthesis = F; 

% Subbands filter lenghts
lf = size(H,1);
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

for i=1:level
    S.AdaptStart(i) = 2^i*S.L(end-i);
end

end

