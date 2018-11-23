clear all;
close all;

% f=audioread('data_test.wav');
% f=f(1:512);
%load noisysig.mat
load leleccum;
f = leleccum;
n = randn(512,1);
order = 3;

t=0:0.01:2;
f=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);


% f = zeros(512,1);
% for i= 1:order
%     f = f + n.^i;
% end;


dwtmode('per','nodisplay');

level = 5;
wname = 'db4';
[C, L] = wavedec(f, level, wname);
D = detcoef(C,L,'cells');
%[cD1, cD2, cD3, cD4, cD5] = detcoef(C,L,[1,2,3,4,5]);
cA = appcoef(C, L, wname);

% %% my_det_coef
% cA2 = C(1:128);
% cD2 = C(129:256);
% cD1 = C(257:512);
% 
% figure;
% subplot(3,1,3)
% plot(dyadup(dyadup(cA2))); title('Level 2 (cA2)'); grid on; axis tight;
% subplot(3,1,2)
% plot(dyadup(dyadup(cA2))); title('Level 2 (cD2)'); grid on; axis tight;
% subplot(3,1,1)
% plot(dyadup(dyadup(cA2))); title('Level 1 (cD1)'); grid on; axis tight;


%%
figure;
subplot(level+2,1,1)
plot(f); title('Original signal');
grid on; axis tight;

WS = 0;
len = length(f);
for i=1:level
    subplot(level+2,1,i+1)
    %plot(D{i});     
    h{i} = dyadup(D{i});
    if i ~= 1   
        for k=1:i-1
            h{i} = dyadup(h{i});
        end
    end
    WS = WS + h{i}(1:len)./2^i;
    plot(h{i});
    title(sprintf('Level %d detail', i));
    grid on; axis tight;
end
subplot(level+2,1,i+2)
plot(cA);
title(sprintf('Level %d apporx', i));
cAA = upsample(cA, 2^level);
WS = WS + cAA(1:len)./2^level;
figure;plot(WS)
