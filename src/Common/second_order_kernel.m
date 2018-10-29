function ker = second_order_kernel(input)
%second_order_kernel
% This function create a random 2nd order kernel for Volterra filtering
% with the most power in first diagonals


if isscalar(input) == 1     % input scalar
    M = input;
    k = rand(M,M)-rand(1);
    k = triu(k);
    flag = 1;
elseif size(input,1) == size(input, 2)  % input symmetric matrix
    k = triu(input);
    M = size(k, 1);
    flag = 1;
else % input vector, build the kernel with it on main diagonal.
    M = numel(input);
    v = input;
    k = zeros(M,M);
    flag = 0;
    
end

if flag == 1

    for i = 0:M-1
        d = diag(ones(M-i,1),i);
        k(d(:,:)==1) = k(d(:,:)==1)./sqrt(i+1);
        
    end
else
    for i = 0:M-1
        d = diag(ones(M-i,1),i);
        shift = floor(i*rand(1));        
        k(d(:,:)==1) = v(shift+1:end-(i-shift))./sqrt(i+1);
    end
end

ker = (k+k') - eye(M).*diag(k);

end

