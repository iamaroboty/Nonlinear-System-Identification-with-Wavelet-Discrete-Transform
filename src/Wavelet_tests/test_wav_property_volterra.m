% test volterra and wavelet optimization using wavelet properties
close all; 
clear all; 


fs = 512;
N = 2*fs;
delta_w = 2/N;
w_freq = (0:N-1)/2; 
w_norm = (0:N-1)*delta_w;
freq = 20;
t = 0:1/fs:1-1/fs;

u = 1*sin(2*pi*freq*t);

figure; 
subplot(3, 1, 1);
plot(t,u);

U = fft(u, N);

subplot(3, 1, 2);
plot(w_freq(1:N/2), abs(U(1:N/2)));
xlabel('Normalized frequency \omega (\pi)'); ylabel('Amplitude');


[low_d,high_d,low_r,high_r] = wfilters('db1');

H = [low_d', high_d'];  % filter matrix analysis
F = [low_r', high_r'];

len = size(low_d,2);

U1 = zeros(len, 2); % linear input 
U2 = zeros(len, 2); % quadratic diagonal
U3 = zeros(len, 2); % x(n)*x(n-1);

er_rr = zeros(len,1);

u1 = zeros(len, 1);
u2 = zeros(len, 1);
u3 = zeros(len, 1);
yn = zeros(size(t,2), 1);

U2_s = zeros(len,2); 
U3_s = zeros(len,2); 

lag =97;

Utap = zeros(len,len,lag);
utap = zeros(lag,1); 

% Analysis and synth
for i=1:size(t,2)    
    u1 = [u(i); u1(1:end-1)];  
    
    u2 = [u(i).*u(i); u2(1:end-1)]; % linear input 
    
    u3 = [u(i).*utap(end); u3(1:end-1)]; % quadratic diagonal term xn*xn
    
    utap = [u(i); utap(1:end-1)]; % off diagnal first x(n)*x(n-1)
    
    % additional len multiplications are performed here
    
    if mod(i,2) == 0         %1 level         
        
        
        
        % method 2 
        
        term1 = u1.*H(:,1); 
        term2 = u1.*H(:,2); 
        
        U_now = [sum(term1); sum(term2)]'; 
        
        U1 = [U_now; U1(1:end-1,:)]; 
        U2 = [u2'*H; U2(1:end-1,:)];
        U3 = [u3'*H; U3(1:end-1,:)];
        
        Utap = cat(3,cat(2,term1 , term2), Utap(:,:,1:end-1)); 
        
        U2_s = [ [sum(term1.*u1); sum(term2.*u1)]' ; U2_s(1:end-1,:)];
        
        loc_lag = floor(lag/2)-1;
        
        if mod(lag,2)==0
        
        
           
        U3_s = [ [sum(Utap(:,1, end-loc_lag).*u1); sum(Utap(:,2, end-loc_lag).*u1)]' ; U3_s(1:end-1,:)]; % shifting property 
        
        else
        
        tmp1 = [Utap(2,1, end-loc_lag-1)  ;Utap(1,1, end-loc_lag)];
        tmp2 = [Utap(2,2, end-loc_lag-1)  ;Utap(1,2, end-loc_lag)];
            
        U3_s = [ [sum(tmp1.*u1); sum(tmp2.*u1)]' ; U3_s(1:end-1,:)];
         
        end
            
      error = (abs(U3)-abs(U3_s)) + (abs(U2)-abs(U2_s));
   
    tmp = (error*F)';
    er_rr = tmp(:,1);
    
    end
    yn(i) = er_rr(1);
    er_rr = [er_rr(2:end); 0];
end

Y = fft(yn, N);

figure; 

plot(t,yn);


       
       
   






