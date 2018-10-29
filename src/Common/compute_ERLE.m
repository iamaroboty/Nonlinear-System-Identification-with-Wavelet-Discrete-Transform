function [ ERLE, indx ] = compute_ERLE( dn , en, n_points )
% plot erle and return fig object 


lendn = max(size(dn)); 

indx = floor(linspace(1000, lendn, n_points)); 
    
ERLE = zeros(n_points,1);
for i = 1:n_points
    
    ERLE(i) = 10*log10((dn(1:indx(i))*dn(1:indx(i))')/(en(1:indx(i))*en(1:indx(i))'));

end


end

