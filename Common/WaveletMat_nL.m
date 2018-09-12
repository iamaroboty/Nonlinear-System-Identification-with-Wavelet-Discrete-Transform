function W = WaveletMat_nL(N,level, varargin)

% WaveletMat_nL      Compute the Wavelet transform matrix for N-levels
%  
% Arguments:
% N                 Signal lenght
% level             levels number
% varargin          wavelet type i.e. 'db1' or 'db4'

if ischar(varargin{1})
    [h,h1,~,~]=wfilters(varargin{1});
else
    h=varargin{1};    h1=varargin{2};
end

W = eye(N);
for i = 1:level
    temp = eye(N);
    temp(1:N/2^(i-1),1:N/2^(i-1)) = WaveletMat_1L(N/2^(i-1),h,h1);
    W = temp*W;
end

