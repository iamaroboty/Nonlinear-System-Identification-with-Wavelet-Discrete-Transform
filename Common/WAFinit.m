function S = WAFinit(w0, mu, level, wtype)

% WSAFinit          Initialize Parameter Structure for Wavelet Transform
%                   Adaptive Filter
%
% Arguments: 
% w0             Coefficients of FIR filter at start (@n=1)
% mu             Step size for LMS algorithm 
% level          Levels of the Wavelet decomposition (dwt)
% wtype          Mother Wavelet filter type

% Assign structure fields
S.coeffs        = w0(:);        % Coefficients of FIR filter 
S.step          = mu;           % Step size 
S.iter          = 0;            % Iteration count
S.AdaptStart    = length(w0);   % Running effects 
M               = length(w0);   % Length of filter
S.beta          = 1/M;          % Forgetting factor
S.levels        = level;        % DWT levels 
S.wtype         = wtype;        % Filter type
S.alpha         = ones(1,M)*1e-6;     % Small positive constant

[h, h1, ~, ~] = wfilters(wtype);
S.W = WaveletMat_nL(M, level, h, h1); % DWT transform matrix
S.analysis = [h', h1'];

end

