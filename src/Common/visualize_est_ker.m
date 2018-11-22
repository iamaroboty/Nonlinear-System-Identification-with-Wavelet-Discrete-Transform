function visualize_est_ker(ker)

lin = ker{1};
quad = ker{2};

if iscell(quad)
    % wavterra
    [M1,~] = size(lin);
    [~,M2] = size(quad);

    k = zeros(M2,M2);
    for i = 0:M2-1
        d = diag(ones(M2-i,1),i);
        k(d(:,:)==1) = quad{i+1};        
    end    
    
else
    % volterra fb
    M2 = (sqrt(8*numel(quad)+1) - 1)  / 2;
    k = zeros(M2,M2);
    
    for i = 0:M2-1
        d = diag(ones(M2-i,1),i);
        
        if i==0
            st = 1;
            en = M2;
        else
            st = en+1;
            en = st+M2-i-1;
        end
            
        k(d(:,:)==1) = quad(st:en);  
    end    
    
k = (k+k') - eye(M2).*diag(k); 

kernel_plot({lin, k})

end