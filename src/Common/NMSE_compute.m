function [ NMSE, indx ] = NMSE_compute( dn , en, n_points )
% plot erle and return fig object 



lendn = max(size(dn)); 

indx = floor(linspace(1000, lendn, n_points)); 
    
NMSE = zeros(n_points,1);
for i = 1:n_points
    
    NMSE(i) = 10*log10((en(1:indx(i))*en(1:indx(i))')/(dn(1:indx(i))*dn(1:indx(i))'));

end


end

