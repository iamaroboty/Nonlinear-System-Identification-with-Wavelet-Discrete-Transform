clear all;
close all;

load wecg;
wtype = 'db3';
level = 3;

f = wecg';

[low_d,high_d,low_r,high_r] = wfilters(wtype);
H = [low_d', high_d']./sqrt(2);  % filter matrix analysis
F = [low_r', high_r']./sqrt(2);  % filter matrix synthesis
[len, ~] = size(H);     % wavelet filter length

x = zeros(length(H),1);
lf = length(H)-1;

delay = (2^level-1)*(length(H)-1)+1;  
fpad = [f, zeros(1,delay)];

for i = 1:level
    delays(i) = 2^i-i;
end

for i = 1:level
    cD{i} = zeros(length(f)+lf*delays(end-i+1),1);
    cA{i} = zeros(length(f)+lf*delays(end-i+1),1);
    bb{i} = zeros(length(F),1);  
end


for n = 1:length(f)+delay
    x = [fpad(n); x(1:end-1)];
    
    tmp = x;
    for i = 1:level        
        xD = H'*tmp;
        cD{i} = [cD{i}(2:end); xD(2)]; 
        cA{i} = [cA{i}(2:end); xD(1)];
        tmp = cA{i}(end:-1:end-lf);        
    end   

    %Synthesis
    for i = level:-1:1
        if i == level
                bb{i} = F*xD + bb{i};            
        else               
                bb{i} = F*[bb{i+1}(1); cD{i}(end-lf*delays(end-i))] + bb{i};
                bb{i+1} = [bb{i+1}(2:end); 0];                        
        end
    end   
    yn(n) = bb{i}(1);
    bb{i} = [bb{i}(2:end); 0];               
    
end

err = max(abs(f-yn(delay:end-1)))