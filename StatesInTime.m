% draw interactive plot in time
figure(11);
clf;
makeplot([], [], out, params);
h = uicontrol('style','slider','units','pixel','position',[20 20 500 20], 'SliderStep', [1e-3 0.1]);
hb = uicontrol('style','pushbutton','units','pixel','position',[540 20 120 20], 'String', 'Rescale Y');
h.Callback = @(hObject, event) makeplot(hObject, event,out, params);
hb.Callback = @(hObject, event) rescalePlot(event);
% addListener(h, 'ContinuousValueChange', @(hObject, event) makeplot(hObject, event,out));

function makeplot(hObject, event,out, params)
if isempty(hObject) 
    ti = 1;
    t = out.t(1);
else
    sval = get(hObject,'Value');
%     t = sval*(out.t(end) - out.t(1)) + out.t(1);
%     ti = find(out.t >= t, 1);
    ti = round(sval*(length(out.t)-1)) + 1;
    t = out.t(ti);
end

% disp(hObject)
% disp(event)

subplot(211);
yyaxis left;cla;hold on;
% Pus = 1 - out.p1_0 - out.p2_0 - out.p3_0;% PU substitute
plot(out.t, out.SL/params.ML,'ro-', out.t, out.LXB/params.ML , '--', 'MarkerSize', 4, 'Linewidth', 2)
% plot(out.t, Pus, 'b-',out.t, out.p1_0, 'r-',out.t, out.p2_0, 'g-',out.t, out.p3_0, 'k-',out.t, out.NR, 'm-');
% plot(out.t, Pus, 'b-', out.t, out.NR, 'm-');

% legend('SL','LXB', 'Pu', 'P1', 'P2', 'P3', 'NR');
% xlim([0 inf])
% ylim([-50, Inf])
plot([t t], [0 1]);
text(t,0.5, num2str(t));
% hold on;
yyaxis right;
plot(out.t, out.Force, 'b-', 'Linewidth', 2)
% legend('SL','LXB', 'Pu', 'P1', 'P2', 'P3', 'NR', 'Force','AutoUpdate','off');
legend('SL','LXB', 'Pu', 'NR', 'Force','AutoUpdate','off');

% plot([simulateTimes;simulateTimes]*1000, repmat([min([fd;out.F']);max([fd;out.F'])], [1 size(simulateTimes, 2)]))
% yyaxis left;


subplot(212);
yl = ylim;
cla;hold on;
% hold on;
% ss = 101; % params.ss;
% s = -0.1:0.002:0.1;%params.s;
ss = params.ss;
% s = params.s + (out.SL(ti) - out.LSE(ti)) - params.LXBpivot;
% s = params.s - (out.SL(ti) - out.LSE(ti)) + params.LXBpivot;
s = params.s' + (-(out.SL(ti) - out.LSE(ti)) + params.LXBpivot)/2;
s = flipud(-s);
% s = params.s + (out.SL(ti) - out.LSE(ti));
% zer = (out.SL(ti) - out.LSE(ti));
zer = 0;
% s_i0 = find(params.s == 0, 1);
% ti = 1;
p1 = out.PU(ti, 1:ss);
p2 = out.PU(ti, ss+1:2*ss);
p3 = out.PU(ti, 2*ss+1:3*ss);
m = max([p1, p2, p3]);
plot(s, p1, '<-b', s, p2, '^-r', s, p3,'>-g', [zer zer], [0 m], '--k');
xlim([s(1), s(end)]);
if ~isequal(yl, [0, 1])
    ylim(yl);
end


end

function rescalePlot(event)
    g = gca;
    ylim('auto');
end