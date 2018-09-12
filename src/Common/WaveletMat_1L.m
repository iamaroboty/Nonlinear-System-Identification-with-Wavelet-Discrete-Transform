function W=WaveletMat_1L(N,varargin)

U=zeros(N/2,N);
V=zeros(N/2,N);
u=zeros(1,N);
v=zeros(1,N);

if ischar(varargin{1})
    [h,h1,~,~]=wfilters(varargin{1});
else
    h=varargin{1};    h1=varargin{2};
end

hrev=wrev(h);
h1rev=wrev(h1);

u(1:length(hrev))=hrev;
v(1:length(h1rev))=h1rev;

if size(U,2)<length(h)           % This may happen when W size is 8 by 8 and
                                 % filter length is more, e.g. bior4.4 (LPF
                                 % of length 9
    U(1,:)=u(1:size(U,2));
    V(1,:)=v(1:size(V,2));
else
    U(1,:)=u;
    V(1,:)=v;
end

for i=2:N/2
      
    U(i,:)=circshift(U(i-1,:)',2)';
    V(i,:)=circshift(V(i-1,:)',2)';
end

W=[U;V];