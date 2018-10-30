function kernel_plot(ker)
% kernel_plot
%
% INPUT ARGUMENT
% ker               Cell array
%
% Plots the first order and second order volterra kernel

assert(iscell(ker), 'Input argument must be a cell array 2x1');

figure('Name', 'Volterra Kernels'); 
subplot(221);   %linear time
plot(ker{1});
ylabel('Linear Kernel (Time)');
xlabel('Sample (N)');
grid on;

subplot(222);   % linear freq
n_points = 1024;
w = linspace(-1,1, n_points);
K1 = fft(ker{1}, n_points);
plot(w, 20*log10(fftshift(abs(K1))));
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Linear Kernel - Power (dB)');
grid on;

sub3 = subplot(223);   % quadratic time
surf(ker{2});
ylabel('Sample (N)');
xlabel('Sample (N)');
view([0 90]);
colormap(sub3, parula);
colorbar
grid on;

sub4 = subplot(224);   % quadratic freq
title('Quadratic Kernel (Frequency)');
K2 = fft2(ker{2}, n_points, n_points);
surf(w, w, 20*log10(fftshift((abs(K2)))), 'LineStyle', 'none');
zlabel('Power (dB)');
ylabel('Normalized Frequency (\times\pi rad/sample)');
xlabel('Normalized Frequency (\times\pi rad/sample)');
view([0 90]);
colormap(sub4, jet);
colorbar
grid on;

end