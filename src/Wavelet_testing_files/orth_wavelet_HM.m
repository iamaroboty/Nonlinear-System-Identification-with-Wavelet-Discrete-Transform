% whitening and orthogonalization property of wavelet trasform 
close all; 
clear all; 

len = 256; 
un = rand(len,1); 

 


ARcoeffs = [5.3217; -9.2948;  7.0933; -2.8152;...
                 2.5805; -2.4230;  0.3747;  2.2628;...
                -0.3028; -1.7444;  1.1053];          % AR(10)
    
% 1/f noise
fv = linspace(0, 1, 20);                                % Normalised Frequencies
a = 1./(1 + fv*2);                                      % Amplitudes Of ‘1/f’
b = firls(42, fv, a);                                   % Filter Numerator Coefficients
%figure(1)
%freqz(b, 1, 2^17)         

 
un = filter(b, 1, un); 


figure;
[cr,lgs] = crosscorr(un,un);
stem(lgs,cr);


% hammerstein pwer series 
order = 5; 

X = zeros(len, order);

for i = 1:order % build the coeff vector (taylor expansion)        
         X(:,i) = un.^i ;                   
end

% plot correlation for power_mat 
  

figure;
    
for row = 1:size(X,2)
    for col = 1:size(X,2)
        [cr,lgs]= crosscorr(X(:,row), X(:,col));
        nm = size(X,2)*(row-1)+col;
        subplot(size(X,2),size(X,2),nm)
         
        stem(lgs,cr,'.');    
        title(sprintf('c_{%d%d}, HM',row,col))
        ylim([0 1])
    end
end




wtype = {'db1'}; 
level = [2]; 

runs = length(wtype)*length(level);
par_comb = combvec(1:length(wtype), 1:length(level));


for i = 1:runs
    fprintf('-------------------------------------------------------------\n');
    fprintf('Run (%d) of (%d)\n', i, runs);           
    
    
    c_level = level(par_comb(2,i));
    c_wtype = wtype{par_comb(1,i)};
    
    
    fprintf('Run hyperpar: wtype = %s, level = %.d \n', c_wtype, c_level);
    
    
  
%     for j = 1:len
% 
%        T = wpdec(X(j,:) , c_level, c_wtype );  % wavelet decomposition 
%        
%        if j==1
%        
%             term_nodes_mat = zeros(size(wpcoef(T,2^c_level-1),2),2.^c_level, len);
%        
%        end
%     
%        for k = 1:2^c_level
%          
%            term_nodes_mat(:,k,j) = wpcoef(T,2^c_level-2+k);
%        end   
%     
%     end
    
    

    T = wpdec2(X , c_level, c_wtype );
    
    
    term_nodes_mat = zeros(size(wpcoef(T,2^c_level-1),1),2.^c_level);
    
    for j = 1:2^c_level
         
       term_nodes_mat(:,j) = wpcoef(T,2^c_level-2+j);
    
    end
    
   X = term_nodes_mat;
    
%    figure;
%           
%    
%    % plot 
%     
%    N = size(X,1);
%    n = (0:N-1)';
% 
%    stem3(n,1:size(X,2),X','filled')
%    ax = gca;
%    ax.YTick = 1:3;
%    view(37.5,30)
%    
%    title(sprintf('wtype: %s ,level: %d',c_wtype, c_level));
    
   figure;
   
    % reduce border effect
   
normval = 0; 
for row = 1:size(X,2)
    for col = 1:size(X,2)
        [cr,lgs] = crosscorr(X(:,row), X(:,col));
        
        nm = size(X,2)*(row-1)+col;
        subplot(size(X,2),size(X,2),nm)
        
        if row ~= col
        normval = normval + norm(cr); 
        end
        
        stem(lgs,cr,'.');    
        title(sprintf('c_{%d%d}, w:%s, lvl:%d',row,col, c_wtype, c_level))
        ylim([0 1])
        
    end
end

fprintf('norm of crosscorr: %.5f \n',normval ); 

    
end


