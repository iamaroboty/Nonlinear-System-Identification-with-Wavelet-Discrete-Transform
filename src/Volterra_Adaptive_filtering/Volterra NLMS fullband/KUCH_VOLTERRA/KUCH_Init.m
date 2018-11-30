function S = KUCH_Init(M, Nc, mu, posdelta)

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
S.mu1          = mu(1);               % Step size 
S.muc        = mu(2);
S.muw         = mu(3);
S.iter          = 0;                % Itertion count 
         % Leaky factor for the leaky LMS algorithm
S.alpha         = 1e-6;             % Small positive constant


S.coeffs{1} = zeros(M(1),1);        % Adaptive filters length
dirac_shifted = zeros(Nc, 1); 
dirac_shifted(posdelta) = 1;
S.coeffs{2} = dirac_shifted; 
S.coeffs{3} = zeros(M(2)-Nc+1,1);   % Nw

S.AdaptStart = max(S.kernel_length);

end
