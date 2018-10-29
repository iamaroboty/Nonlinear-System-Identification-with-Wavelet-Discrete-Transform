%clear all;

f=audioread('data_test.wav');
f=f(1:512);
%load noisysig.mat
% load leleccum;
% f = leleccum;

%dwtmode('per','nodisplay');

level = 5;
wname = 'db1';
[C, L] = wavedec(f, level, wname);
D = detcoef(C,L,'cells');
%[cD1, cD2, cD3, cD4, cD5] = detcoef(C,L,[1,2,3,4,5]);
cA5 = appcoef(C, L, wname);
%% my_det_coef
cA2 = C(1:128);
cD2 = C(129:256);
cD1 = C(257:512);

figure;
subplot(3,1,3)
plot(dyadup(dyadup(cA2))); title('Level 2 (cA2)'); grid on; axis tight;
subplot(3,1,2)
plot(dyadup(dyadup(cD2))); title('Level 2 (cD2)'); grid on; axis tight;
subplot(3,1,1)
plot(dyadup(cD1)); title('Level 1 (cD1)'); grid on; axis tight;


%%
figure;
subplot(level+1,1,1)
plot(f); title('Original signal');
grid on; axis tight;

for i=1:level
    subplot(level+1,1,i+1)
    %plot(D{i}); 
    h{i} = dyadup(D{i});
    if i ~= 1   
        for k=1:i-1
            h{i} = dyadup(h{i});
        end
    end
    plot(h{i});
    title(sprintf('Level %d detail', i));
    grid on; axis tight;
end
