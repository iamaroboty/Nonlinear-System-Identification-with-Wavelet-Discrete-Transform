function [Hi,H_af] = create_petraglia_structure(H,level)

Hi = create_multilevel_bank(H,level);

indx = 1; 
indx2 = 1 ; 

for i= 1:size(Hi,2)                          
    H_af(:,indx) = conv(Hi(:,i), Hi(:,i));
    indx = indx +1; 

    if i+1 <= size(Hi,2)
        H_af(:,indx) = conv(Hi(:,i), Hi(:,i+1));
        indx = indx +1; 
    end             
end

end

