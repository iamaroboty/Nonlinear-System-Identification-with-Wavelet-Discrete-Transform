function [Hi] = create_multilevel_bank(H,level)
%CREATE_MULTILEVEL_BANK 

if level == 1
    Hi = H;
else
    indx = 1;
    for i = 1:size(H,2)
        for j=1:level-1
            up{i,j} = upsample(H(:,i), 2^(j)); 
        end
    end
    
    %outer product
    H_tmp = H; 
    for i=1:size(up,2)
        H_tmp = outer_conv(H_tmp, up(:,i));
    end
    Hi = H_tmp;
    
%     %calculate useless z in H
%     nz = 0;
%     for i = 0:(level-1)
%         nz = 2^i + nz;
%     end
%     nz = nz - level;
%   
%     Hi = H_tmp(1:end-nz,:);
        
end

