%% DISCRETE WAVELET DECOMPOSITION TESTING
%
% Testing Signal
clear all; close all

d = 256;
t=0:0.001:1;
un=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);  
un = un(1:d);

level = 1;


%% Generazione dei coefficienti del filtro
df = dtfilters('dtf1');

%% Decomposizone con funzione matlab
wt = dddtree('cplxdt',un,level,df{1},df{2});

reconstructed_matlab = idddtree(wt);


%% Decomposizione brutale teorica
cmod = 'full';


fdf = df{1};
df = df{2};

low_d_first_1 = fdf{1}(:,1); 
high_d_first_1 = fdf{1}(:,2); 

low_d_first_2 = fdf{2}(:,1); 
high_d_first_2 = fdf{2}(:,2);

low_d_1 = df{1}(:,1); 
high_d_1 = df{1}(:,2); 

low_d_2 = df{2}(:,1); 
high_d_2 = df{2}(:,2);


% Banco di analisi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level 1
d_cA_1=conv(un,low_d_first_1,cmod);
d_cA_1=downsample(d_cA_1,2,1);

d_cD_1=conv(un,high_d_first_1,cmod);
d_cD_1=downsample(d_cD_1,2,1);
[McA_1,McD_1] = dwt(un,low_d_first_1,high_d_first_1);

% Level 2
d_cA2_1=conv(d_cA_1,low_d_1,cmod);
d_cA2_1=downsample(d_cA2_1,2,1);

d_cD2_1=conv(d_cA_1,high_d_1,cmod);
d_cD2_1=downsample(d_cD2_1,2,1);
[McA2_1,McD2_1] = dwt(McA_1,low_d_1,high_d_1);



%% second tree

d_cA_2=conv(un,low_d_first_2,cmod);
d_cA_2=downsample(d_cA_2,2,1);

d_cD_2=conv(un,high_d_first_2,cmod);
d_cD_2=downsample(d_cD_2,2,1);
[McA_2,McD_2] = dwt(un,low_d_first_2,high_d_first_2);

% Level 2
d_cA2_2=conv(d_cA_2,low_d_2,cmod);
d_cA2_2=downsample(d_cA2_2,2,1);

d_cD2_2=conv(d_cA_2,high_d_2,cmod);
d_cD2_2=downsample(d_cD2_2,2,1);
[McA2_2,McD2_2] = dwt(McA_2,low_d_2,high_d_2);

stem(wt.cfs{1,1}(:,:,1)); 
hold on; 
stem(McD_1); 



%% Test Real time
len = size(df{1}{1},1);
M= d; 
L = [M; zeros(level,1)];
L = [M; zeros(level,1)];
for i= 1:level
    L = [floor((L(1)+len-1)/2); L(1:end-1)];
end
S.L = [L(1); L]';

fdf = df{1};
df = df{2};

% H1{1} = fdf{1};
% H2{1} = fdf{2};
% 
% F1{1} = flipud(fdf{1});
% F2{1} = flipud(fdf{2});
% 
% H1{2} = df{1};
% H2{2} = df{2};
% 
% F1{2} = flipud(df{1});
% F2{2} = flipud(df{2});


for i= 1:level
    
    % filters
    if level == 1
    H1{i} =  fdf{1};
    H2{i} = fdf{2};
    F1{1} = flipud(fdf{1});
    F2{1} = flipud(fdf{2});
    
    
    end
    
    H1{i} = df{1};
    H2{i} = df{2};
    F1{i} = flipud(df{1});
    F2{i} = flipud(df{2});
    
    
    U1.cD{i} = zeros(L(end-i),1);    
    U1.cA{i} = zeros(L(end-i),1);   
    U2.cD{i} = zeros(L(end-i),1);    
    U2.cA{i} = zeros(L(end-i),1);    
    Y1.cD{i} = zeros(L(end-i),1);
    Y1.cA{i} = zeros(L(end-i),1);
    Y2.cD{i} = zeros(L(end-i),1);
    Y2.cA{i} = zeros(L(end-i),1);
    
    eD1{i} = zeros(L(end-i),1);      % Error signa, transformed domain
    eDr1{i} = zeros(len,1);          % Error signal, time domain
    
    
    
    eD2{i} = zeros(L(end-i),1);      % Error signa, transformed domain
    eDr2{i} = zeros(len,1);          % Error signal, time domain
    
    delays(i) = 2^i-1;              % Level delay for synthesis
    w1{i} = zeros(L(end-i),1);       % Subband adaptive filter coefficient, initialize to zeros    
    w1{i}(1) = 1;
    
    w2{i} = zeros(L(end-i),1);
    w2{i}(1) = 1; 
    
end 
w1{i} = zeros(L(end-i),2);           % Last level has 2 columns, cD and cA
w1{i}(1,1)=1; 
w1{i}(1,2)=1; 
eD1{i} = zeros(1,2);                 % Last level has 2 columns, cD and cA
U1.tmp = zeros(len,1);
Y1.tmp = zeros(len,1);
U1.Z = zeros(2,1);
Y1.Z = zeros(2,1);

w2{i} = zeros(L(end-i),2);           % Last level has 2 columns, cD and cA
w2{i}(1,1)=1; 
w2{i}(1,2)=1; 

eD2{i} = zeros(1,2);                 % Last level has 2 columns, cD and cA
U2.tmp = zeros(len,1);
Y2.tmp = zeros(len,1);
U2.Z = zeros(2,1);
Y2.Z = zeros(2,1);

%pwr = w;
%beta = 1./L(2:end-1);

u = zeros(len,1);                 % Tapped-delay line of input signal (Analysis FB)  
y = zeros(len,1);                 % Tapped-delay line of desired response (Analysis FB)


ITER = length(un);
en = zeros(1,ITER);               % Initialize error sequence to zero
iter=0;


for n = 1:ITER    
    u = [un(n); u(1:end-1)];        % Input signal vector contains [u(n),u(n-1),...,u(n-M+1)]'
      

    % Analysis Bank
    U1.tmp = u;
    
    
    U2.tmp = u;
    
    
    for i = 1:level
        if mod(n,2^i) == 0
    
            U1.Z = H1{i}'*U1.tmp;
            U1.cD{i} = [U1.Z(2); U1.cD{i}(1:end-1)]; 
            U1.cA{i} = [U1.Z(1); U1.cA{i}(1:end-1)];
            U1.tmp = U1.cA{i}(1:len);
            
            U2.Z = H2{i}'*U2.tmp;
            U2.cD{i} = [U2.Z(2); U2.cD{i}(1:end-1)]; 
            U2.cA{i} = [U2.Z(1); U2.cA{i}(1:end-1)];
            U2.tmp = U2.cA{i}(1:len);
            
                   
            
            if i == level
                eD1{i} = sum(([U1.cA{i}, U1.cD{i}]).*w1{i});
                
                eD2{i} = sum(([U2.cA{i}, U2.cD{i}]).*w2{i});
                

            else
                eD1{i} = [eD1{i}(2:end); U1.cD{i}'*w1{i}]; 
                
                
                eD2{i} = [eD2{i}(2:end); U2.cD{i}'*w2{i}]; 
                
            
            end           
            iter = iter + 1;                
        end
    end    


    % Synthesis Bank
    for i = level:-1:1     
        if i == level
            if mod(n,2^i) == 0
                eDr1{i} = F1{i}*eD1{i}' + eDr1{i};
                eDr2{i} = F2{i}*eD2{i}' + eDr2{i};
                
            end
        else
            if mod(n,2^i) == 0                
                eDr1{i} = F1{i}*[eDr1{i+1}(1); eD1{i}(end-(len-1)*delays(end-i))] + eDr1{i};
                eDr1{i+1} = [eDr1{i+1}(2:end); 0];
                
                eDr2{i} = F2{i}*[eDr2{i+1}(1); eD2{i}(end-(len-1)*delays(end-i))] + eDr2{i};
                eDr2{i+1} = [eDr2{i+1}(2:end); 0];
            end            
        end
    end   
    en(n) = 0.5*(eDr1{i}(1) + eDr2{i}(1));
    eDr1{i} = [eDr1{i}(2:end); 0];       
    eDr2{i} = [eDr2{i}(2:end); 0];  
end

en = en(1:ITER);

