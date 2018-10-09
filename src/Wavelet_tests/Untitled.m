% Calculate the boundary filters for an M-channel orthogonal case.
%
% Cormac Herley, Columbia University
%
% For details see "Boundary Filters for Finite Length Signals
% and Time-varying Filter Banks," C. Herley, submitted 1993.
%
% Expects H0, H1, ... H_{M-1} as row vectors, produces, a square unitary
% matrix Q, which contains boundary filters at each end.
%
M = 2;
block = [];
ka = floor(length(H0-1)/M); na=length(H0);
p = ka+1;
for i=1:M
text=sprintf('block = [block; H%g]',i-1);
eval(text);
end
%
D=block*block';
block=sqrt(1/D(1,1))*block;
% Calculate the truncated (rank-deficient) matrix
Ht=zeros(M*p,M*(p+ka-1));
for i=1:ka
for j=1:M
mat=block(1:M,(i-1)*M+1:M*i);
end
repeat=[zeros(p,(i-1)) eye(p) zeros(p,ka-i)];
Ht=Ht+ kron(repeat,mat);
end
% Do Gram-Schmidt; use a random vector as input
% Of course this just produces any old set of
% boundary filters. Should optimise these.
da=0;
Q = Ht; sz = M*(p+ka+da-1); nsd = da+ka-1;

% Get the left boundary filters first; keep going until
% the nullspace is spanned.
dl=0; cnd=1.0;
while (cnd-1.0)^2 < 1.0e-5
tr = rand(sz-na);
left = [tr(1,:) zeros(1,na)];
lt = left*(eye(sz) - Q'*Q);
nrm=norm(lt,2);
cnd = cond([lt/norm(lt,2); Q]);
if (cnd-1.0)^2 < 1.0e-5
Q = [lt/norm(lt,2); Q];
dl=dl+1;
end
end
% Remaining vectors are in the right nullspace.
dr=0; cnd=1.0;
while (cnd-1.0)^2 < 1.0e-5
tr = rand(sz-na);
right=[zeros(1,na) tr(2,:)];
rt = right*(eye(sz) - Q'*Q);
cnd = cond([Q; rt/norm(rt,2)]);
if (cnd-1.0)^2 < 1.0e-5
Q = [Q; rt/norm(rt,2)];
dr=dr+1;
end
end