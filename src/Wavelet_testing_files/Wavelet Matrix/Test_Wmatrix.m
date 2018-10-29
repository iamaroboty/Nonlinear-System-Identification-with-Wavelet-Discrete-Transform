%Wavelet discrete decomposition, normal and matrix form
clear all; close all

filtername='db2'; 
levels = 2; 
len = 128;

t=0:0.0001:1;
f=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
f = f(1:len*2);
%f = [zeros(len-1,1); 1; zeros(len,1)]'; 

N=length(f); 

[low_d,high_d,low_r,high_r] = wfilters(filtername);
W=WaveletMatrix_nL(N,levels,filtername) ;

delta32_dwt = (W*f')';
f_r = (W'*delta32_dwt')';

figure; plot([f',f_r'])

% for dval=1:levels 
%     D{dval} = delta32_dwt((N/(2^(dval)))+1:(N/(2^(dval-1)))); 
%     if (dval == levels) % only last level for cA
%         A{dval} = delta32_dwt(1:(N/(2^(dval)))); 
%     end 
% end

% dwtmode('zpd');
[C, L] = wavedec(f, levels, filtername);
% cD = detcoef(C,L,'cells');
% cA = appcoef(C, L, filtername);

%figure; plot([C', delta32_dwt']);
