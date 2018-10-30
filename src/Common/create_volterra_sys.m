function [Sys_obj] = create_volterra_sys(order, lengths, gains, name )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if order ~= size(lengths,2)
    
    error('order and lengths should have same dimension!');
    
end
Sys_obj.order = order;
Sys_obj.name = name; 

for i = 1:order    
    Sys_obj.M(i) = lengths(i); % append lengths of kernels 
    
    if order > 2
        
        error("not supported order >2");
        
    end
    
    if i ==1
        
    
    randFIR = randIR(lengths(i)) ;    
        
    Sys_obj.Responses{i} = (gains(i).*randFIR);
    
    elseif i ==2
        
    tmp = triu(ones(lengths(i))); 
        
    for j = 0:lengths(i)-1
        
        
       
            
        randFIR = randIR(lengths(i)-j); 
        
        d = diag(ones(lengths(i)-j,1),j);
        tmp(d(:,:)==1) = tmp(d(:,:)==1).*randFIR;

    end
    
     
                           
    Sys_obj.Responses{i} =   (tmp+tmp')' ; % if matrix is square this is symmetric 
    
    
    end
    
end

end

