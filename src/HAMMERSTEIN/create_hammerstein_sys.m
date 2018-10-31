function [System] = create_hammerstein_sys(M, Gains)
%create_kernel
% Possible modes: "delta", "randomdiag", "random", "lowpass", "simulated"
% THIS CODE NEEDS TO ME OPTIMIZED, BUT THIS WORKS.



for i = 1:size(M,2)
    
   System.filters{i} = Gains(i).*rand(M(i), 1);  
    
end

System.M = M; 


end


