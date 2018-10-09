%Wavelet discrete decomposition, normal and matrix form
clear all; close all

filtername='db4'; 
levels = 8; 
len = 256;

t=0:0.0001:1;
f=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
f = f(1:len*2);
%f = [zeros(len-1,1); 1; zeros(len,1)]'; 

N=length(f); 

[low_d,high_d,low_r,high_r] = wfilters(filtername);
[C, L] = wavedec(f, levels, filtername);

W1 = WMatrix1(N, filtername);   % first level transform matrix

%first level decomposition
V1 = (W1*f')';  
[C1, L1] = wavedec(f, 1, filtername);
cA = V1(1:N/2);
sumerr1 = sum((V1 - C1).^2)

N2 = N/2;
assert(mod(N2,2)==0, 'N2 must be even')

W2 = WMatrix1(N2, filtername);
V2 = (W2*cA')';

V = [V2 V1(N2+1:end)];
sumerr = sum((V-C).^2)

T = eye(N);
T(1:N2,1:N2) = W2;
W = T*W1;   % 2 level matrix transform

VV = (W*f')';