
 clear all; close all
% Segnale di prova%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:0.001:2;
f=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
%f = [sin(2*pi*10*t) sin(2*pi*100*t)];
%f=sin(2*pi*t*5000)+sin(2*pi*t*150)+sin(2*pi*t*1050);
% f=randn(1,1000);
% f=f/max(abs(f));
%f=audioread('data_test.wav');
f=f(1:256);
%f= f(1:end-2);
d=length(f);
m=1:d/2;
% audiowrite('data_test.wav',f,44100);
%figure; plot(f);
%%
% Generazione dei coefficienti del filtro %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[low_d,high_d,low_r,high_r] = wfilters('db2');

% Decomposizone con funzione matlab
[C, L] = wavedec(f, 2, 'db2');

%Conv mode
cmod = 'full';

% Banco di analisi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level 1
d_cA=conv(f,low_d,cmod);
d_cA=downsample(d_cA,2,1);

d_cD=conv(f,high_d,cmod);
d_cD=downsample(d_cD,2,1);
[McA,McD] = dwt(f,low_d,high_d);

% Level 2
d_cA2=conv(d_cA,low_d,cmod);
d_cA2=downsample(d_cA2,2,1);

d_cD2=conv(d_cA,high_d,cmod);
d_cD2=downsample(d_cD2,2,1);
[McA2,McD2] = dwt(McA,low_d,high_d);

% Banco di sintesi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level 2 
r_cA2=upsample(d_cA2,2);
r_cA2=conv(r_cA2,low_r,cmod);

r_cD2=upsample(d_cD2,2);
r_cD2=conv(r_cD2,high_r,cmod);

r_l2=r_cA2+r_cD2;
r_l2=r_l2(1:end-1);
McR2 = idwt(McA2,McD2,low_r,high_r);

% Level 1 
r_cA=upsample(r_l2,2);
r_cA=conv(r_cA,low_r,cmod);

r_cD=upsample(d_cD,2);
r_cD=conv(r_cD,high_r,cmod);

out=r_cA+r_cD;
out=out(1:end-1);
McR = idwt(McR2,McD,low_r,high_r);

figure;
plot(f);
hold on;
plot(McR,'r');
hold on
plot(out,'b');

%figure;


% plot(cd1);hold on;plot(cd11);
% plot(cd2);hold on;plot(cd22);figure;
% plot(cd3);hold on;plot(cd33);figure;
% % a1=filter(1,low_d,f);lowpass2=downsample(a1,2);lowpass2=lowpass2*2;
% % b1=filter(1,high_d,f(2:end));highpass2=downsample(b1,2);highpass2=highpass2/2;
% 
% a3=f(2*m-1).*v(1) + f(2*m).*v(2);
% d3=f(2*m-1).*w(1) + f(2*m).*w(2);
% 
% au2=upsample(ad1,2);
% a2=filter(1,HI_R,au2);
% 
% 
% bu2=upsample(bd1,2);
% b2=filter(1,LO_R,bu2);
% 
% 
% out=a2+b2;
% % a2=downsample(a1,2);
% % d2=downsample(b1,2);
% % 
% % a2=upsample(a1,2);
% % b2=upsample(b1,2);
% 
% 
% 
% 
% 
% % out=out/max(abs(out));
% 
% plot(out);
% %a2=a2*2;
% % a1=f(2*m-1).*v(1) + f(2*m).*v(2);
% % plot(a1); hold on; plot(a2);
% % figure;
% % d1=f(2*m-1).*w(1) + f(2*m).*w(2);
% % plot(d1);hold on; plot(d2);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %plots
% % figure(1)
% % subplot(3,2,[1 2])
% % plot(f)
% % title('Original Signal') %plot of the original signal
% % ylabel('Amplitude')
% % xlabel('Time')
% axis tight
% 
% subplot(3,2,3)
% plot(a1)
% title('First Trend') %plot of first trend after haar wavelet transform
% ylabel('Amplitude')
% xlabel('Time')
% axis tight
% 
% subplot(3,2,4)
% plot(d1)
% title('First Fluctuation') %plot of first fluctuation after haar wavelet transform
% ylabel('Amplitude')
% xlabel('Time')
% axis tight
% 
% subplot(3,2,[5 6])
% plot([a1 d1])
% title('Transformed Signal') %plot of the final signal obtained after the transformation
% xlabel('Amplitude')
% ylabel('Time')
% axis tight

%end