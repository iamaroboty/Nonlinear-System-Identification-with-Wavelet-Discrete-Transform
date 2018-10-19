%%
clear all;
close all;

level = 2;
wtype = 'db16';
DFTpoint = 1024;
delta_w = 2*pi/DFTpoint;
w = (0:DFTpoint-1)*delta_w;
w2 = w(1:DFTpoint/2);

[low_d,high_d,low_r,high_r] = wfilters(wtype);
H = [low_d', high_d'];     % Analysis filters
F = [low_r', high_r'];     % Synthesis filters

[Hi, H_af] = create_petraglia_structure(H, level);

%% Check orthogonality of the bank and QMF condition
m = fft(H(:,1), DFTpoint);
mt = fft(H(:,2), DFTpoint);

figure;
subplot(221); plot(w2/pi,20*log10(abs(m(1:512)))+eps); 
title(sprintf('Transfer modulus of %s', wtype))
axis([0 1 -100 10])
xlabel('\omega (\pi)'); ylabel('Power (dB)');
subplot(222); plot(w2/pi,20*log10(abs(mt(1:512)))+eps); 
title(sprintf('Transfer modulus of QMF %s', wtype))
axis([0 1 -100 10])
xlabel('\omega (\pi)'); ylabel('Power (dB)');
subplot(2,2, [3 4]); plot(w/pi,abs(m).^2 + abs(mt).^2); 
axis([0 2 1.999 2.001])
title(sprintf('Check QMF condition for %s and QMF %s', wtype, wtype)) 
xlabel(sprintf(' abs(fft(%s))^2 + abs(fft(qmf(%s))^2 = 1', wtype, wtype))

% Check for orthonormality. 
df = [H(:,1),H(:,2)]; 
id = df'*df;

if (id - eye(2)) < 1e-5
    fprintf('Filters %s are orthogonal\n', wtype);
else
    fprintf('Filterbank is not orthogonal\n');
end

%% 
HHi = fft(Hi, DFTpoint);
HH_af = fft(H_af, DFTpoint);

% Plot full level bank
figure; subplot(3,1, [1 2]); hold on;
for i = 1:size(Hi')
    plot(w2/pi, 20*log10(abs(HHi(1:DFTpoint/2, i))+eps), 'LineWidth', 2,...
        'Displayname', ['H' num2str(i-1)]);
end
title(sprintf('Complete (%d level) filterbank, %s', level, wtype));
ylabel('Power (dB)');
legend('show'); grid on;
axis([0 1 -100 inf]); box on;

% Compute and plot distortion function
THi = sum(abs(HHi').^2);
subplot(3,1,3); plot(w2/pi,10*log10(THi(1:DFTpoint/2)));
xlabel('Frequency,\omega (\pi)'); ylabel('Gain (dB)');
% axis([0 1 -0.01 0.01]);

% Plot petraglia structure
figure; hold on;
for i = 1:size(H_af')
    if mod(i,2) == 0
        plot(w2/pi, 20*log10(abs(HH_af(1:DFTpoint/2, i))+eps), 'LineWidth', 2, ...
            'Displayname', ['H' num2str(floor((i-1)/2)) 'H' num2str(i/2)]);
    else
        plot(w2/pi, 20*log10(abs(HH_af(1:DFTpoint/2, i))+eps), 'LineWidth', 2, ...
            'Displayname', ['H' num2str(floor((i-1)/2)) 'H' num2str(floor((i-1)/2))]);
    end
end
title(sprintf('Petraglia structure'));
xlabel('\omega (\pi)'); ylabel('Power (dB)');
legend('show'); grid on;
axis([0 1 -100 inf])

fprintf(' Sum power of H: \n');
disp(sum(abs(fft(H)).^2));
fprintf(' Sum power of Hi (fullbank): \n');
disp(sum(abs(fft(Hi)).^2));
fprintf(' Sum power of Petragliabank: \n');
disp(sum(abs(fft(H_af)).^2));

% Sum of petraglia for synthesis
Ti = sum(abs(HHi).^2);
Taf = sum(abs(HH_af).^2);
%Powers in petraglia subbbands
Tsum = [Taf(1)+Taf(2) , Taf(2)+Taf(3)+Taf(4), Taf(4)+Taf(5)+Taf(6), Taf(6)+Taf(7)]

