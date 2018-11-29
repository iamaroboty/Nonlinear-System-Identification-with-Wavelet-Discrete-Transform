function [out] = visualize_est_ker(ker)

lin = ker{1};
quad = ker{2};

if iscell(quad)
    % wavterra
    [M1,~] = size(lin);
    [~,M2] = size(quad);
    [kerlen, ~] = size(quad{1});

    k = zeros(kerlen,kerlen);
    for i = 0:M2-1
        d = diag(ones(kerlen-i,1),i);
        k(d(:,:)==1) = quad{i+1};        
    end    
    
    M2 = kerlen;
    
else
    % volterra fb
    M2 = (sqrt(8*numel(quad)+1) - 1)  / 2;
    k = zeros(M2,M2);
    
    for i = 1:M2
        d = [zeros(1,i-1),ones(1,M2-i+1)];
        
        if i==1
            st = 1;
            en = M2;
        else
            st = en+1;
            en = st+M2-i;
        end
            
        k(d==1,i) = quad(st:en);  
    end    
    
end

k = (k+k') - eye(M2).*diag(k); 

kernel_plot({lin, k})

if nargout > 0
    out = k;
end

end