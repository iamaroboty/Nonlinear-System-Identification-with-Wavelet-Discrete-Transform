%
clear all;
close all;

d = 256;        %Total signal length
M = 4;
wtype = 'db1';
t=0:0.001:10;
un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
% un = 1:d;
un = un(1:d)';

[low_d,high_d,low_r,high_r] = wfilters(wtype);
H = [low_d', high_d'];  % filter matrix analysis
F = [low_r', high_r'];  % filter matrix synthesis
[len, ~] = size(H);     % wavelet filter length

index = [];                     % Indexing for correct second order input
for i=0:M-1
    index = cat(2,index, M*i+1:M*(i+1)-i);
end

for i = 1:M
    w{i} = zeros(M, 1);
    WcA{i} = zeros(M/2);
    WcD{i} = zeros(M/2);
end
    
utmp = zeros(M^2-M+1,1);
cA = zeros(floor((M^2-M+1)/2));
cD = zeros(floor((M^2-M+1)/2));
cA1 = zeros(floor((M^2-M+1)/2));
cD1 = zeros(floor((M^2-M+1)/2));
u = zeros(M,1);

for n = 1:length(un)
    u = [un(n); u(1:end-1)]; 
                 
    utmp = [un(n)*u; utmp(1:end-M)];    % Temporal delay line
    u2 = utmp(index);  
    
    for i = 1:M
        w{i} = [un(n)*u(i); w{i}(1:end-1)];
    end
    
    if mod(n,2) == 0
        U2 = H'*u2(1:len);
        cA = [U2(1), cA(1:end-1)];
        cD = [U2(2), cD(1:end-1)];
        
        for i = 1:M
            W{i} = H'*w{i}(1:len);
            WcA{i} = [W{i}(1), WcA{i}(1:end-1)];
            WcD{i} = [W{i}(2), WcD{i}(1:end-1)];
        end
    end
end
        
    
    