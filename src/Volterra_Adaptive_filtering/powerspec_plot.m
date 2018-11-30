function powerspec_plot(AR, iter)
% powerspec_plot

% AR coeffs
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
ARcoeffs(6).a = [1; -0.9];  
ARcoeffs(7).a = [1; 0.9]; 


DFTpoint = 4096;
delta_w = 2*pi/DFTpoint;
w = (0:DFTpoint-1)*delta_w;
w2 = w(1:DFTpoint/2);

wn = randn(iter,1);   

un = filter(1,ARcoeffs(AR).a,wn);            % IIR filtering
HzSig = fft(ARcoeffs(AR).a, DFTpoint);
HzSig = HzSig';                 % Magnitude response of AR coloring filter
HzSig = 1./(abs(HzSig).^2); 
HzSig = HzSig/var(un);          % Normalized to unit variance

figure('Name', 'PowerSpectrum');
plot(w2/pi, 10*log10(HzSig(1:DFTpoint/2)+eps), 'DisplayName', ['AR = ', num2str(AR)], ...
            'LineWidth', 1.5);
ylabel('Power (dB)');
xlabel('\omega (\pi)');
grid on; legend('show');

end

