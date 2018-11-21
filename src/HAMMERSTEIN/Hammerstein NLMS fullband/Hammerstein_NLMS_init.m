function S = Hammerstein_NLMS_init(order, M, mu, leaks, alpha)

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
S.filters_lengths  = M;              % Kernel dimensions
S.step          = mu;             % Step saize 
S.alpha         = alpha;           % Small positive constant
S.leaks         = leaks; 
S.iter          = 0; 
S.AdaptStart    = M;
S.order         = order; 

for j=1:M,                    % DCT-transform matrix
   S.T(1,j)=1/sqrt(M);
   for i=2:M,
     S.T(i,j)=sqrt(2/M)*cos(pi*(i-1)*(2*j-1)/2/M);
   end
end	

end
