% Test delay line
clear all;
close all;

un = 1:200;
K = [2 4];
L = K(2)*(K(2)+1)/2;

for i = 1:K(2)
    D{i} = zeros(K(2)-i+1, 1);
end

u = zeros(K(2),1);

 
for n = 1:length(un)
    u = [un(n); u(1:end-1)]; 
    
    for i = 1:K(2)
        D{i} = [un(n)*u(i); D{i}(1:end-1)];
    end

    u2 = cat(1, D{1:end});

end
