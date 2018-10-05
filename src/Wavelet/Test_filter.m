clear all; close all;
a = [ 1 1 1]; 
b = [1 2 3];
x = 1:10;

MA = filter(b,1,x);
AR = filter(1,a,x);
ARMA = filter(b,a,x)

%% Conv with Toepliz matrix
r = [b zeros(1,length(x)-length(b))];
c = [b(1) zeros(1,length(x)-1)];
xConv = toeplitz(c,r);
xt = x*xConv;
xc = conv(x,b);
diff = max(abs(xt-xc(1:length(x))))

%% Conv with for
xpad = [zeros(1,length(b)-1), x];
xfor = zeros(1, length(x));
z = zeros(1,length(b));
for i = 1:length(x)
%     xfor(i) = flip(b)*xpad(i:length(b)+i-1)'
    xfor(i) = b(1) * x(i) + z(1);
   for m = 2:length(b)
      z(m - 1) = b(m) * x(i) + z(m);
   end
end
diff_for = max(abs(MA-xfor(1:length(x))))


%% AR-MA (Direct Form 2 II - Transpose)
bn = b / a(1);
an = a / a(1);
z = zeros(1,length(a));
Y = zeros(size(x));
for m = 1:length(Y)
   Y(m) = bn(1) * x(m) + z(1);
   for i = 2:length(a)
      z(i - 1) = bn(i) * x(m) + z(i) - an(i) * Y(m);
   end
end
z = z(1:length(a) - 1);

diff_my = max(abs(ARMA - Y))

