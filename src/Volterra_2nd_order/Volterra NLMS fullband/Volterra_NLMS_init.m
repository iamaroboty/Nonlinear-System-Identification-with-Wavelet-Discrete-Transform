function S = Volterra_NLMS_init(M, mu, leak)

% Volterra_NLMS_init          Initialize Parameter Structure Volterra
%                               Fullband filtering
% Arguments: 
% M              Kernel lengths
% mu             Step size for LMS algorithm 

if nargin < 3                    % Set default to conventional NLMS
	leak = 0;
end

% Assign structure fields
S.kernel_length  = M;               % Kernel dimensions
S.step          = mu;               % Step size 
S.iter          = 0;                % Itertion count 
S.leakage       = leak;          % Leaky factor for the leaky LMS algorithm
S.alpha         = 1e-6;             % Small positive constant

M1 = M(1);
M2 = M(2)*(M(2)+1)/2;       % Combination coefficient
S.length = [M1, M2];        % Adaptive filters length

S.AdaptStart = max(S.kernel_length);

end
