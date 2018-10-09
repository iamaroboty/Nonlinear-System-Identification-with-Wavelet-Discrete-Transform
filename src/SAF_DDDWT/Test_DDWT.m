%% DOUBLE DENSITY DWT TESTING

% Testing Signal
clear all; close all

d = 256;        %Total signal length
t=0:0.001:10;
f=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
f = f(1:d)';
% f=f(1:256);
%f = [1; -10; 324; 48; -483; 4; 7; -5532; 34; 8889; -57; 54];
%d=length(f);

% d = 512;
% f = load('h1.dat');         % Unknown system (select h1 or h2)
% f = f(1:d);                 % Truncate to length M

% f = zeros(d,1);
% f(3) = 1;

% f = chirp((0:0.001:2),0,1,250);
% f = f(1:d);

type = 'ddt';
filters = 'filters1';
level = 1;

%% MATLAB FUNCTION
wt = dddtree(type,f,level,filters);
% wt.cfs{2} = zeros(1,512);
xrec = idddtree(wt);
% plot(t(1:d),xrec,'linewidth',2)
% set(gca,'xtick',[0 0.3 0.72 1]); set(gca,'xgrid','on');
err = max(abs(f-xrec))

%% Real time
[H, F] = dtfilters(filters);
[len, ~]= size(H);
delay = (2^level-1)*(len-1); 

L = [d; zeros(level,1)];
for i= 1:level
    L = [d/(2^i); L(1:end-1)];
    delays(i) = 2^i-1;
end
L = [L(1); L]'; 

for i= 1:level
    cA{i} = zeros(L(end-i),1);
    cD{i} = zeros(L(end-i),2);
    bb{i} = zeros(len,1);    
end

x = zeros(len,1);
fpad = [f; zeros(delay,1)];
yn = zeros(d+delay, 1);
xD = zeros(3,1);

for n = 1:d+delay
    x = [fpad(n); x(1:end-1)];
    
    tmp = x;
    for i = 1:level
        if mod(n,2^i) == 0
          xD = H'*tmp;
          cA{i} = [cA{i}(2:end); xD(1)];
          cD{i}(:,1) = [cD{i}(2:end,1); xD(2)]; 
          cD{i}(:,2) = [cD{i}(2:end,2); xD(3)]; 
          tmp = cA{i}(end:-1:end-(len-1));
        end
    end
    
    % Synthesis
    for i = level:-1:1
        if i == level
            if mod(n,2^i) == 0
                bb{i} = F*xD + bb{i};
            end
        else
            if mod(n,2^i) == 0                
                bb{i} = F*[bb{i+1}(1); cD{i}(end-(len-1)*delays(end-i),1); cD{i}(end-(len-1)*delays(end-i),2);] + bb{i};
                bb{i+1} = [bb{i+1}(2:end); 0];
            end            
        end
    end   
    yn(n) = bb{i}(1);
    bb{i} = [bb{i}(2:end); 0];                   
end

err_DDDWT = max(abs(f-yn(1+delay:end))) 