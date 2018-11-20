N = 11;
n = (0:N-1)';

a = 0.4;
b = 0.7;
c = 0.999;

xabc = [a.^n b.^n c.^n];
figure;
stem3(n,1:3,xabc','filled')
ax = gca;
ax.YTick = 1:3;
view(37.5,30)

waitforbuttonpress

%%
[cr,lgs] = xcorr(xabc,'coeff');

for row = 1:3
    for col = 1:3
        nm = 3*(row-1)+col;
        subplot(3,3,nm)
        stem(lgs,cr(:,nm),'.')
        title(sprintf('c_{%d%d}',row,col))
        ylim([0 1])
    end
end