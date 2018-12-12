%% FIGURE MODIFY SCRIPT
% Input MATLAB figure as .fig

clear all;

%Insert here figure name
%or just open and set the figure to be modified as current figure
uiopen('SOAF_comp.fig',1)

% General values
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',0.8);
set(0,'defaultLineMarkerSize',5);

% Getting variables to workspace, double click on them and GUI pops up!!
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% objTypes = get(dataObjs, 'Type');  %type of low-level graphics object

%For a beamer slide: width=5.04 in, length=3.78 in
%Values in inch
width = 20;
heigth = 6;

set(h,...    
    'Units', 'centimeters', ...
    'Position', [1 1 width heigth],...
    'PaperPositionMode', 'auto')

set(gca,...
    'Units', 'normalized',...
    'Position', [.15 .2 .75 .7],...
    'FontUnits', 'points',...
    'FontWeight', 'normal',...
    'FontSize', 8,...
    'FontName', 'Serif')

%% custom values, depends on the figue
ylim([-50 10])
title('Input signal: colored with AR(1)')

%% PRINT TO .EPS
print('SOAF_comp.eps','-depsc2')
