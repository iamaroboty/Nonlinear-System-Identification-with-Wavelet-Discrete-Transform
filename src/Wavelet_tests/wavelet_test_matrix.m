%% DISCRETE WAVELET DECOMPOSITION TESTING
%
% Testing Signal
clear all; close all

d = 256;
t=0:0.001:1;
f=20*(t.^2).*(1-t).^4.*cos(12*t.*pi)+sin(2*pi*t*5000)+sin(2*pi*t*150);
f = f(1:d)';
% f=f(1:256);
%f = [1; -10; 324; 48; -483; 4; 7; -5532; 34; 8889; -57; 54];
%d=length(f);

% d = 512;
% f = load('h1.dat');         % Unknown system (select h1 or h2)
% f = f(1:d);                 % Truncate to length M

wtype = 'db1';
level = 2;

%% Generazione dei coefficienti del filtro
[low_d,high_d,low_r,high_r] = wfilters(wtype);
W = WaveletMat_nL(d, level, low_d, high_d); % DWT transform matrix
H = [low_d', high_d'];  % filter matrix analysis
F = [low_r', high_r'];  % filter matrix synthesis
[len, ~] = size(H);     % wavelet filter length

%% Decomposizone con funzione matlab
%dwtmode('asymh')
[C, L] = wavedec(f, level, wtype);
f_rec = waverec(C, L, wtype);

%% Decomposizione con Matrice W
Z = W*f;
Zr = W'*Z;

for dval=1:level 
    D{dval} = Z((d/(2^(dval)))+1:(d/(2^(dval-1)))); 
    if (dval == level) % only last level for cA
        A{dval} = Z(1:(d/(2^(dval)))); 
    end 
end


%% Decomposizione brutale teorica
cmod = 'full';
% Banco di analisi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Level 1
d_cA=conv(f,low_d,cmod);
d_cA=downsample(d_cA,2,1);

d_cD=conv(f,high_d,cmod);
d_cD=downsample(d_cD,2,1);
[McA,McD] = dwt(f,low_d,high_d);

% Level 2
d_cA2=conv(d_cA,low_d,cmod);
d_cA2=downsample(d_cA2,2,1);

d_cD2=conv(d_cA,high_d,cmod);
d_cD2=downsample(d_cD2,2,1);
[McA2,McD2] = dwt(McA,low_d,high_d);

C2 = [McA2; McD2; McD];

% % Banco di sintesi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Level 2 
r_cA2=upsample(d_cA2,2);
r_cA2=conv(r_cA2,low_r,cmod);

r_cD2=upsample(d_cD2,2);
r_cD2=conv(r_cD2,high_r,cmod);

r_l2=r_cA2+r_cD2;
r_l2=r_l2(1:end-1);
 McR2 = idwt(McA2,McD2,low_r,high_r);

% Level 1 
r_cA=upsample(r_l2,2);
r_cA=conv(r_cA,low_r,cmod);

r_cD=upsample(d_cD,2);
r_cD=conv(r_cD,high_r,cmod);

out=r_cA+r_cD;
out=out(1:end-1);
 McR = idwt(McR2,McD,low_r,high_r);

%% Old testing
% for i=1:d    
%     a = [f(i); a(1:end-1)];  
%     A = [a, A(:,1:end-1)];
%     if mod(i,2)== 0         %1 level         
%         U1 = (H'*A)';           
%         U2 = U(1:end-2,:);
%         U = [U1', U2']';   % U = [cA1, cD1] not decimated
%         
%         ca = [U(1,1), ca(1:end-1)]';
%         cA = [ca, cA(:,1:end-1)];
%         if mod(i,4) == 0        %2 level
%             V1 = (H'*cA)';
%             V2 = V(1:end-2,:);
%             V = [V1', V2']';    % V = [cA2, cD2] not decimated
%         end
%     end
% end

a = zeros(len,1);
U = zeros(d/2,2);
V = zeros(d/4,2);
Z = zeros(d/2,1);
fr = zeros(d,1);
% Analysis and synth
for i=1:d    
    a = [f(i); a(1:end-1)];  
    if mod(i,2)== 0         %1 level         
        U = [a'*H; U(1:end-1,:)]; %U(:,1) = cA, U(:,2) = cD %both flipped upside down
    if mod(i,4) == 0
        V = [U(1:len)*H; V(1:end-1,:)]; %V(:,1) = cA2, V(:,2) = cD2 %both flipped upside down
        
        Z = [flip(V(1,:)*F)'; Z(1:end-2,:)];
        fr = [flip([Z(2), U(2,2)]*F)'; fr(1:end-2,:)];
    elseif i>=4
        fr = [flip([Z(1), U(2,2)]*F)'; fr(1:end-2,:)];
    end
    end
    if i == d
        fr = [flip([Z(1), U(1,2)]*F)'; fr(1:end-2,:)];
    end
end

% %Synthesis
z = zeros(d/4,1);
for i=0:d/4-1
    z = [flip(V(end-i,:)*F)'; z(1:end-1, :)]; %z = McR2
end

f_r = zeros(d/2,1);
for i=0:d/2-1
    f_r = [flip([z(end-i), U(end-i,2)]*F)'; f_r(1:end-1, :)];
end

lev1 = zeros(len, 2);
lev2 = zeros(1, 2);
lev2_r = zeros(len,1);
x = zeros(len, 1);
fr = zeros(d, 1);
% Analysis and synth
for i=1:d    
    x = [f(i); x(1:end-1)];  
    if mod(i,2)== 0         %1 level         
        lev1 = [x'*H; lev1(1:end-1,:)]; %U(:,1) = cA, U(:,2) = cD %both flipped upside down
    if mod(i,4) == 0
        lev2 = lev1(1:len)*H;
        
        lev2_r =  flip(lev2(1,:)*F)';
        lev1_r = flip([lev2_r(2), lev1(2,2)]*F)';
    elseif i>=4
        lev1_r = flip([lev2_r(1), lev1(2,2)]*F)';
    end
    end
    
end

%% DWT Real Time
lev1 = zeros(len, 2);
lev2 = zeros(1, 2);
lev2_rr = zeros(len,1);
lev1_rr = zeros(len,1);
x = zeros(len, 1);
fr = zeros(d, 1);
% Analysis and synth
for i=1:d    
    x = [f(i); x(1:end-1)];  
    if mod(i,2)== 0         %1 level         
        lev1 = [x'*H; lev1(1:end-1,:)]; %U(:,1) = cA, U(:,2) = cD %both flipped upside down
    if mod(i,4) == 0
        lev2 = lev1(1:len)*H;
        
        lev2_rr =  (lev2(1,:)*F)';
    end
    lev1_rr = ([lev2_rr(1) , lev1(2,2)]*F)';
    lev2_rr = [lev2_rr(2:end); 0];    
    end
    yn(i) = lev1_rr(1);
    lev1_rr = [lev1_rr(2:end); 0];
end

%% DTW Test Real time
yn = zeros(d,1);
lev1 = zeros(d/2, 2);
lev2 = zeros(d/4, 2);
lev2_rr = zeros(d/2,1);
lev1_rr = zeros(len,1);
x = zeros(len, 1);
fr = zeros(d, 1);
% Analysis and synth
for i=1:d    
    x = [f(i); x(1:end-1)];  
    if mod(i,2)== 0         %1 level         
        lev1 = [x'*H; lev1(1:end-1,:)]; %U(:,1) = cA, U(:,2) = cD %both flipped upside down
    if mod(i,4) == 0
        lev2 = [lev1(1:len)*H; lev2(1:end-1,:)];
        
        lev2_rr =  (lev2(1,:)*F)' + lev2_rr(1,:)';
    end
    lev1_rr = ([lev2_rr(1) , lev1(2,2)]*F)';
    lev2_rr = [lev2_rr(2:end); 0];    
    end
    yn(i) = lev1_rr(1);
    lev1_rr = [lev1_rr(2:end); 0];
end

% %% DTW Test Others Wavelets %% Got some problems
% yn = zeros(d,1);
% lev1 = zeros(d/2, 2);
% lev2 = zeros(d/4, 2);
% lev2_rr = zeros(d/2,1);
% lev1_rr = zeros(len,1);
% x = zeros(len, 1);
% fr = zeros(d, 1);
% % Analysis and synth
% for i=1:d    
%     x = [f(i); x(1:end-1)];  
%     if mod(i,2)== 0         %1 level         
%         lev1 = [x'*H; lev1(1:end-1,:)]; %U(:,1) = cA, U(:,2) = cD %both flipped upside down
%     if mod(i,4) == 0
%         lev2 = [lev1(1+1:len+1)*H; lev2(1:end-1,:)];
%         
%         lev2_rr =  (reshape(lev2([1+2,2+2],:),1,4)*F)' + lev2_rr(1,:)';
%     end
%     lev1_rr = ([lev2_rr([1,2]) ; lev1([3+1,4+1],2)]'*F)';
%     lev2_rr = [lev2_rr(2:end); 0];    
%     end
%     yn(i) = lev1_rr(1);
%     lev1_rr = [lev1_rr(2:end); 0];
% end
    
% figure;
% plot(f);
% hold on;
% plot(McR,'r');
% figure;
% plot(wrev(f_r));
% hold on
% plot(out,'b');