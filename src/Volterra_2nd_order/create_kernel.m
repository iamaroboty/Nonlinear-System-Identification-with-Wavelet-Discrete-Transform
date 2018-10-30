function [ker1, ker2] = create_kernel(M1, M2, mode, param)
%create_kernel
% Possible modes: "delta", "randomdiag", "random", "lowpass", "simulated"
% THIS CODE NEEDS TO ME OPTIMIZED, BUT THIS WORKS.

rng('default');

fprintf('2nd order Volterra kernel. Type: "%s" \n', mode);
% Kronecker delta in deltapos
switch mode
    case 'delta'
        deltapos = param{1};
        ker1 = zeros(M1,1);
        ker1(deltapos(1)) = 1;
        ker2 = zeros(M2,M2);
        ker2(deltapos(2),deltapos(2)) = 1;
        fprintf('Delta position: \n');
        disp(deltapos);
    
        
    case 'randomdiag'
        ker1 = rand(M1,1) - rand(1);
        shift = 0;      
        ker2 = diag(rand(M2-shift,1)- rand(1) , shift); 

        N = param{2}; %diagonals number, beyond the main one
        for i = 1:N
            d = diag(ones(M2-i,1),i);
            ker2(d(:,:)==1) = rand(M2-i, 1) - rand(1);
        end    

        d = eye(M2); ker2(d(:,:)==1) = rand(M2,1)- rand(1) ;     % instert principal diagonal
        fprintf('Number of diagonals: \n');
        disp(N);
    
        
    case 'random'
        ker1 = rand(M1,1)-rand(1);
        ker2 = second_order_kernel(M2);
       
        
    case 'lowpass'
        normfreq = param{3};
        if numel(normfreq) == 1
            normfreq(2) = normfreq(1);
        end
        samples = [M1 M2]/2-1;
        b1 = normfreq(1)*sinc(normfreq(1)*(-samples(1)-1:samples(1)));
        b2 = normfreq(2)*sinc(normfreq(2)*(-samples(2)-1:samples(2)));
        ker1 = b1;
        ker2 = second_order_kernel(b2);
        fprintf('Normalized cutoff frequency: \n');
        disp(normfreq);
    
        
    case 'simulated'
        h1h2 = param{4};
        b1 = load(h1h2(1));
        b1 = b1(1:M1);
        ker1 = b1;

        b2 = load(h1h2(2));
        b2 = b2(1:M2);
        ker2 = second_order_kernel(b2);
        
        disp(h1h2);
        
        
    otherwise
        modes = ["delta", "randomdiag", "random", "lowpass", "simulated"];
        error('mode must be one of those chars: %s', sprintf('%s ', modes));
        
end


