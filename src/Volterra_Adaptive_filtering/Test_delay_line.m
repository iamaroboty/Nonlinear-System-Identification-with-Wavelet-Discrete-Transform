% Test delay line
clear all;
close all;

un = 1:200;
K = [2 5];
L = K(2)*(K(2)+1)/2;

for i = 1:K(2)
    D{i} = zeros(K(2)-i+1, 1);
end

M = K(2);
index = [];
for i=0:M-1
    index = cat(2,index, M*i+1:M*(i+1)-i);
end

u = zeros(M,1);
uu2 = zeros(M*M-M+1, 1);

 
for n = 1:length(un)
    u = [un(n); u(1:end-1)]; 
    
    for i = 1:K(2)
        D{i} = [un(n)*u(i); D{i}(1:end-1)];
    end

    u2 = cat(1, D{1:end});
    uu2 = [un(n).*u; uu2(1:end-K(2))];    
    x2 = uu2(index);

end
