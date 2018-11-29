function [est] = Volterra_test(un, ker)

ker_quad = visualize_est_ker(ker);
ker_lin = ker{1};
[M1,~] = size(ker_lin);
[~,M2] = size(ker{2});

est = fastVMcell(un,{ker_lin, ker_quad}, [M1, M2] );
est = sum(est,1);

end