% reenact_exp_script.m
% Function: f(x) = 1 - A * exp(-(x/s).^2)

clearvars; close all; clc

% Parameters
A0 = 0.08;   % depth
s0 = 3;      % width
x = linspace(-10,10,2000);
% f = @(x,A,s) 1 + A .* exp(-(x./s)) + s*0;
f = @(x,A, s) min(1, (x<0).*exp(-A.*x) + (x>=0).*(1 - exp(-A.*(x+1))))+ s*0;
f = @(x,A, s) min(1, 1 - 1./(1 + exp(-A.*(1-x))))+ s*0;
f = @(x,A, s) min(1, (x>=0).*(1 - exp(-A.*x)) + (x<0).*(1./(1 + exp(A.*x))))+ s*0;
% vectorized lambdas
u = @(x,s) min( max( (1 - x)./s, 0 ), 1 );    % clamped parameter in [0,1]
smoothstep = @(t) 3.*t.^2 - 2.*t.^3;                  % C1 smoothstep on [0,1]

% main function: f(x,k,delta)
% - k > 0 controls steepness (larger k = steeper drop near 0)
% - delta > 0 controls how quickly the tail is forced to 1 around x=1
f = @(x,A,s) 1 + ( exp(-A.*x) - exp(-A) ) .* smoothstep( u(x,s) );
f = @(x,A,s) 1 + (A - 1).*(1 - x).^s;

f = @(x,A,s) (x <= 1).* (1 + (A - 1).*(1 - x).^s) + (x > 1);




% Figure + axes
hFig = figure('Name','Exponential reenact','NumberTitle','off',...
    'Units','normalized','Position',[0.2 0.2 0.6 0.55]);
hAx = axes('Parent',hFig,'Position',[0.08 0.22 0.9 0.72]); 
hold(hAx,'on'); grid(hAx,'on');
xlabel(hAx,'x'); ylabel(hAx,'f(x)');
title(hAx,'f(x) = 1 - A * exp(-(x/s).^2)');

% Plot curve
hPlot = plot(hAx, x, f(x,A0,s0),'LineWidth',2);

% Markers and text
xPts = [-5 0 5];
hPts = plot(hAx, xPts, f(xPts,A0,s0),'o','MarkerFaceColor','w','LineWidth',1.5);
hTxt = gobjects(1,numel(xPts));
for k=1:numel(xPts)
    hTxt(k) = text(xPts(k), f(xPts(k),A0,s0)+0.005,...
        sprintf('x=%g: %.4f',xPts(k), f(xPts(k),A0,s0)),...
        'HorizontalAlignment','center','FontSize',10,'Parent',hAx);
end

% Info text
hInfo = uicontrol('Style','text','Units','normalized','Position',[0.02 0.02 0.96 0.14],...
    'FontSize',10,'HorizontalAlignment','left','BackgroundColor',get(hFig,'Color'));

% Sliders
uicontrol('Style','text','Units','normalized','Position',[0.05 0.16 0.12 0.03],'String','A (depth)');
hA = uicontrol('Style','slider','Units','normalized','Position',[0.17 0.16 0.36 0.05],...
    'Min',0.01,'Max',2,'Value',A0,'SliderStep',[0.001 0.01]);

uicontrol('Style','text','Units','normalized','Position',[0.55 0.16 0.12 0.03],'String','s (width)');
hs = uicontrol('Style','slider','Units','normalized','Position',[0.67 0.16 0.36 0.05],...
    'Min',0.1,'Max',20,'Value',s0,'SliderStep',[0.01 0.1]);

% Store handles in the figure (so callbacks can access them)
handles = struct('x',x,'f',f,'xPts',xPts,'hPlot',hPlot,'hPts',hPts,'hTxt',hTxt,...
                 'hInfo',hInfo,'hAx',hAx,'hA',hA,'hs',hs);
guidata(hFig,handles);

% Assign callbacks
hA.Callback = @(src,evt)updatePlot(hFig);
hs.Callback = @(src,evt)updatePlot(hFig);

% Initial update
updatePlot(hFig);

% --------- Callback function (not nested) ----------
function updatePlot(hFig)
    handles = guidata(hFig);
    A = handles.hA.Value;
    s = handles.hs.Value;
    y = handles.f(handles.x,A,s);
    set(handles.hPlot,'YData',y);

    yPts = handles.f(handles.xPts,A,s);
    set(handles.hPts,'YData',yPts);

    for ii=1:numel(handles.xPts)
        handles.hTxt(ii).Position(2) = yPts(ii)+0.005;
        handles.hTxt(ii).String = sprintf('x=%g: %.4f',handles.xPts(ii),yPts(ii));
    end

    handles.hInfo.String = sprintf('A=%.4f, s=%.2f | f(-5)=%.4f, f(0)=%.4f, f(5)=%.4f',...
        A,s,yPts(1),yPts(2),yPts(3));
    % ylim(handles.hAx,[min(y)-0.02 1.01]);
    drawnow;
end
