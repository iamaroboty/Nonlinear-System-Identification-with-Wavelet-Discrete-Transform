function [Y, z] = myFilter(b, a, X, z)

na  = length(a);
nb  = length(b);

if na > nb
    b = [b, zeros(1,na-nb)];
elseif nb > na
    a = [a, zeros(1,nb-na)];
end
n = max(na,nb);
z(n) = 0;      % Creates zeros if input z is omitted, z is the delay history of the filter

% Coefficient normalization
b = b / a(1);  
a = a / a(1);

Y    = zeros(size(X));
for m = 1:length(Y)
   Y(m) = b(1) * X(m) + z(1);
   for i = 2:n
      z(i - 1) = b(i) * X(m) + z(i) - a(i) * Y(m);
   end
end
z = z(1:n - 1);
end