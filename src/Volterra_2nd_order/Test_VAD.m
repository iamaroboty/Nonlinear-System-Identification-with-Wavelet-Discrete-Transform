% Testing VAD
clear all;
close all;

addpath(genpath('../Common')); 

filename = 'timit_m.mat';  

load(filename,'un');

un = un./max(un);
fs = 8000;

frame = 80;
Eth = 40;

u = zeros(frame,1);                

ITER = length(un);
VAD = zeros(1,ITER);

Efr = [];


for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'   

    if mod(n,80) == 0   % Each frame, no overlap (10 ms)
        Efr = cat(2, Efr, 10*log10(sum(u.^2)));
        Emin = min(Efr);
        Eth = Eth+(Emin);
    end
    
    if n > 80
        if (10*log10(sum(u.^2)) - Emin) >= Eth
            VAD(n) = 1;
        end
    end                          
end

t = 1:ITER;
figure; plot(t, un);
hold on; plot(t(VAD == 1), un(VAD == 1), 'r');

