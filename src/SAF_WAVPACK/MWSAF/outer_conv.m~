function [ H ] = outer_conv( vector1, vector2 )

H = zeros(size(vector2{1},1)+size(vector1,1)-1, size(vector2,1)*size(vector1,2)); 

indx = 1;
for i=1:size(vector1,2)
   
    for j=1:size(vector2,1)
       
        H(:, indx) = conv(vector1(:,i), vector2{j});
   
        indx = indx +1; 
        
        
              
    end
    
    
end


end

