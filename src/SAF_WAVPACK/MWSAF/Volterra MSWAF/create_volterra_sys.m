function [Sys_obj] = create_volterra_sys(order, lengths, name )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if order ~= size(lengths,1)
    
    error('order and lengths should have same dimension!');
    
end
Sys_obj.order = order;
Sys_obj.name = name; 

for i = 1:order    
    Sys_obj.M = lengths(i); % append lengths of kernels 
    Sys_obj.Responses{i} = rand(lengths(i),1);
end

end

