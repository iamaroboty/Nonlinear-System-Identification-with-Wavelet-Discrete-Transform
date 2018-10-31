function S = KUCH_NLMS_init(M, Nc, mu)

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
         % Leaky factor for the leaky LMS algorithm
S.alpha         = 1e-6;             % Small positive constant


S.coeffs{1} = zeros(M(1),1);        % Adaptive filters length
dirac_shifted = zeros(Nc, 1); 
dirac_shifted(1) = 0.5;
S.coeffs{2} = dirac_shifted; 
S.coeffs{3} = zeros(M(2)-Nc+1,1);

S.AdaptStart = max(S.kernel_length);

end
