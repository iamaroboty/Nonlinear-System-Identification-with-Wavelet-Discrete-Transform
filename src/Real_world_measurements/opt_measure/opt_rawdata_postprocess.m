clear all
close all

noise = -130; %dBfs

[dn, ~] = audioread('behr_ref_color.wav');
dn = dn(:,1);
[raw, ~] = audioread('opt_colored.wav');
raw = raw(:,1);

fs = 44100;
est_delay = 104;    %ms
est_delay_sample = round(fs*104e-3);

figure; 
subplot(311); plot(dn); 
subplot(312); plot(raw);

un = raw(1:length(dn));
delay = finddelay(un,dn);
delay = abs(delay);

d = delay;

un = un(d:end);
dn = dn(1:end-d+1);
assert(length(dn)==length(un));

subplot(313); plot(un);

save('opt_color.mat', 'un');
save('opt_ref_color.mat', 'dn');


[P, Q] = rat(8000/fs);
dnn = resample(dn, P, Q);
unn = resample(un, P, Q);