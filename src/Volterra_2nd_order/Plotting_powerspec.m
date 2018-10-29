clear all;
close all;

%% AR coeffs
ARcoeffs(1).a = 1;                                   % White noise
ARcoeffs(2).a = [1; -0.1; -0.8];                     % AR(2), c.f. Figure 6.14(a)
ARcoeffs(3).a = [1; -1.6;  0.81];                    % AR(2), c.f. Figure 4.6
ARcoeffs(4).a = [5.3217; -9.2948;  7.0933; -2.8152;...
                 2.5805; -2.4230;  0.3747;  2.2628;...
                -0.3028; -1.7444;  1.1053];          % AR(10), c.f. Figure 6.17
ARcoeffs(5).a = [1.0000; -1.3193;  0.8610; -0.4541;...
                -0.6191;  0.9800; -0.6530;  0.5424;...
                 0.3694; -0.5102;  0.3673; -0.4017;...
                 0.0649; -0.1499;  0.1212;  0.1978]; % AR(15), c.f. Figure 3.3(a)
ARcoeffs(6).a = [1; 0.9];  

%%
DFTpoint = 4096;
delta_w = 2*pi/DFTpoint;
w = (0:DFTpoint-1)*delta_w;
w2 = w(1:DFTpoint/2);

randn('state',sum(100*clock));  % Reset random generator to diffrent state
wn = randn(4*80000,1);          % white Gaussian random signal, unit variance
ar = 1:6;


%% Spectrum of colored signal
figure;
for i=1:6
    un = filter(1,ARcoeffs(i).a,wn);            % IIR filtering

    HzSig = fft(ARcoeffs(i).a, DFTpoint);
    HzSig = HzSig';                 % Magnitude response of AR coloring filter
    HzSig = 1./(abs(HzSig).^2); 
    HzSig = HzSig/var(un);          % Normalized to unit variance

    subplot(6,1,i);
    plot(w2/pi, 10*log10(HzSig(1:DFTpoint/2)+eps), 'DisplayName', ['AR = ', num2str(i)], ...
            'LineWidth', 1.5);
    ylabel('Power (dB)');
    grid on;
    legend('show');
    
end
xlabel('\omega (\pi)');


%%
M1 = 256; % length of first order volterra kernel
M2 = 32; % length of second order volterra kernel

b1 = load('h1.dat');
b1 = b1(1:M1);
ker1 = b1;

b2 = load('h2.dat');
b2 = b2(1:M2);
ker2 = second_order_kernel(b2);

% ker1 = rand(M1,1)-rand(1);
% ker2 = second_order_kernel(M2);

NL_system.M = [M1, M2];
NL_system.Responses = {ker1, ker2};

kernel_plot(NL_system.Responses);

fs = DFTpoint/2;
freq_w = pi/4;
freq = fs*freq_w / (2*pi);
dt = 1/fs;
un = 1*sin(2*pi*freq*(0:dt:5-dt));

max_len= max(NL_system.M); 
dn = fastVMcell(un, NL_system.Responses , NL_system.M);
dn = sum(dn,1);

% Adjusting starting index of signals

dn = dn(1:length(un));      
dn = dn(max_len+1:end);
un = un(max_len+1:end);

% Normalization of u(n) and d(n)

un = un/std(dn);
dn = dn/std(dn);

UN = fft(un, DFTpoint);
DN = fft(dn, DFTpoint);


figure;
subplot(211);
plot(w2/pi, 20*log10(abs(UN(1:DFTpoint/2)+eps)), 'DisplayName', 'Power of un');
ylabel('Power (dB)');
subplot(212);
plot(w2/pi, 20*log10(abs(DN(1:DFTpoint/2)+eps)), 'DisplayName', 'Power of dn');
ylabel('Power (dB)');
xlabel('\omega (\pi)');
legend('show');


%%
% figure;
% subplot(211);
% plot(w2/pi, 20*log10(abs(UN(1:DFTpoint/2))), 'LineWidth', 1.5, ...
%     'DisplayName', 'Power Spectrum of un')
% ylabel('Power (dB)');
% grid on; 
% legend('show');
% 
% subplot(212);
% plot(w2/pi, 20*log10(abs(DN(1:DFTpoint/2))), 'LineWidth', 1.5, ...
%     'DisplayName', 'Power Spectrum of dn')
% xlabel('\omega (\pi)');
% ylabel('Power (dB)');
% grid on; 
% legend('show');
% axis([0 1 -100 inf]); 