% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;
clear all;

% Defaults for this blog post
width = 5.04;     % Width in inches
height = 3.78;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

f = @(x) x.^2;
g = @(x) 5*sin(x)+5;
dmn = -pi:0.001:pi;
xeq = dmn(abs(f(dmn) - g(dmn)) < 0.002);

%%
% % The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

% Now we repeat the first example but do not need to include anything
% special beyond manually specifying the tick marks.
figure(1); clf;
plot(dmn,f(dmn),'b-',dmn,g(dmn),'r--',xeq,f(xeq),'g*');
xlim([-pi pi]);
legend('f(x)', 'g(x)', 'f(x)=g(x)', 'Location', 'SouthEast');
xlabel('x');
title('Automatic Example Figure');
% set(gca,'XTick',-2:0); %<- Still need to manually specific tick marks
set(gca,'YTick',0:10); %<- Still need to manually specific tick marks

%%
% figure('Units','inches',...
%         'Position',[0 0 width height],...
%         'PaperPositionMode','auto');
% 
% t = 0:.1:4*pi;
% y = sin(t);
% plot(t,y)
% xlabel('Time(s)')
% ylabel('y(t)')
% title('Sin function')
% legend('y=sin(t)')    
%     
% axis([0 t(end) -1.5 1.5])
% set(gca,...
% 'Units','normalized',...
% 'YTick',-1.5:.5:1.5,...
% 'XTick',0:t(end)/4:t(end),...
% 'Position',[.15 .2 .75 .7],...
% 'FontUnits','points',...
% 'FontWeight','normal',...
% 'FontSize',9,...
% 'FontName','Times')
 
    