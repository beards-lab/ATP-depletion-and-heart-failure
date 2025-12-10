% draw interactive plot in time
t_init = 2.77;
t_init = 0;
t_init = out.t(round(end/2));

figure(10);clf;
% params = getParams(params);
plot(out.t, out.p1_0, '-', out.t, out.p2_0, '-', out.t, out.PuATP, '-',out.t, out.PuR, '-', out.t, out.SR , '-', out.t, out.SRD, '-', LineWidth=1.5)
legend('P1','P2','PuATP','PuR', 'SRT', 'SRD')
fig = figure(11);
clf;
% makeplot([], [], out, params);
h = uicontrol('style','slider','units','pixel','position',[65 20 500 20], 'SliderStep', [1e-3 0.1]);
lbl = uicontrol('Style','edit','units','pixel','Position',[5 20 55 20]); % <----- #1
hb = uicontrol('style','pushbutton','units','pixel','position',[570 20 120 20], 'String', 'Rescale Y');

% set init value
uisliderSetVal(h, out.t, t_init);

h.Callback = @(hObject, event) makeplot(hObject, event,out, out.params, lbl);
hb.Callback = @(hObject, event) rescalePlot(event);
lbl.Callback = @(hObject, event) setTime(hObject, event,out, out.params, h);
% addListener(h, 'ContinuousValueChange', @(hObject, event) makeplot(hObject, event,out));

makeplot(h, [], out, out.params, lbl);

function setTime(uiedit, ~,out, params, uislider)
    val = str2num(uiedit.String);

    if isempty(val)
        return;
    end

    uisliderSetVal(uislider, out.t, val);
    makeplot(uislider, [],out, params, []);
end

function uisliderSetVal(uislider, t_arr, t)
    ti = find(t_arr>t, 1);
    uislider.Value  = (ti - 1)/(length(t_arr)-1);
end

function makeplot(hObject, ~,out, params, lbl)
%%
if ~exist("hObject","var") || isempty(hObject) 
    % ti = 1;
    ti = find(out.t>2.8, 1);
    ti = length(out.t);
else
    sval = get(hObject,'Value');
%     t = sval*(out.t(end) - out.t(1)) + out.t(1);
%     ti = find(out.t >= t, 1);
    ti = round(sval*(length(out.t)-1)) + 1;
end
t = out.t(ti);

if exist("lbl","var") 
    lbl.String = t;
end

% disp(hObject)
% disp(event)
subplot(2,3,1:2);
yyaxis left;cla;hold on;
% Pus = 1 - out.p1_0 - out.p2_0 - out.p3_0;% PU substitute
plot(out.t, out.SL/params.ML,'ro-', out.t, out.LXB/params.ML , '-', 'MarkerSize', 4, 'Linewidth', 2)
% plot(out.t, Pus, 'b-',out.t, out.p1_0, 'r-',out.t, out.p2_0, 'g-',out.t, out.p3_0, 'k-',out.t, out.NR, 'm-');
% plot(out.t, Pus, 'b-', out.t, out.NR, 'm-');
% plot(out.t, out.NR, 'k-');

% legend('SL','LXB', 'Pu', 'P1', 'P2', 'P3', 'NR');
% xlim([0 inf])
% ylim([-50, Inf])
plot([t t], [0 1]);
printtext = sprintf('%0.6fs \n%0.2f um\n%0.2f kPa\n%0.2f s^{-1}', t, out.SL(ti), out.Force(ti), out.RTD(ti));
text(t,0.5, printtext);
% hold on;
yyaxis right;
plot(out.t, out.Force, 'b-', 'Linewidth', 2)
% legend('SL','LXB', 'Pu', 'P1', 'P2', 'P3', 'NR', 'Force','AutoUpdate','off');
legend('SL','LXB', 'Force','AutoUpdate','off');

% plot([simulateTimes;simulateTimes]*1000, repmat([min([fd;out.F']);max([fd;out.F'])], [1 size(simulateTimes, 2)]))
% yyaxis left;
subplot(2, 3, [3 6]);cla;
plotStateFluxes(out, t);

subplot(2,3,4:5);
yl = ylim;
% cla;hold on;

% hold on;
% ss = 101; % params.ss;
% s = -0.1:0.002:0.1;%params.s;
ss = params.ss;
% s = params.s + (out.SL(ti) - out.LSE(ti)) - params.LXBpivot;
% s = params.s - (out.SL(ti) - out.LSE(ti)) + params.LXBpivot;
% s = params.s' + (-(out.SL(ti) - out.LSE(ti)) + params.LXBpivot)/2;
if params.LegacyStrainFlipping
    s = params.s' + (-(out.SL(ti) - out.LSE(ti)) + out.LXBPivot(ti))/2;
    s = flipud(-s);
else
    s = params.s - (-(out.SL(ti) - out.LSE(ti)) + out.LXBPivot(ti))/2;
end
% s = params.s + (out.SL(ti) - out.LSE(ti));
% zer = (out.SL(ti) - out.LSE(ti));
zer = 0;
% s_i0 = find(params.s == 0, 1);
% ti = 1;
dS = params.dS;
p1 = out.PU(ti, 1:ss)*dS;
p2 = out.PU(ti, ss+1:2*ss)*dS;
if params.NumberOfStates == 3
    p3 = out.PU(ti, 2*ss+1:3*ss)*dS;
else
    p3 = nan(size(p2));
end


% p1_0 = dS*sum(p1);% p1_1 = dS*sum(s.*p1);
% p2_0 = dS*sum(p2); 
% p2_1 = params.dS*sum(s'.*p2);
% p3_0 = dS*sum(p3); 

% if params.UseP31Shift
%     p3_1 = params.dS*sum((s'+params.dr).*p3);
%     % p3_1_stroke = dS*sum((s+params.dr).*p3.*(s+params.dr >= 0));
%     % p3_1_overstroke = dS*sum((s+params.dr).*p3.*(s+params.dr < 0));
% 
% else
%     p3_1 = params.dS*sum(s'.*p3);
% end

m = max([p1, p2, p3]);
% plot(s, p1, '<-b', s, p2, '^-r', s+params.dr, p2, '--r', s+params.dr, p3,'>-g', [zer zer], [0 m], '--k', LineWidth=1.5);
plot(s, p1, '<-b', s, p2, '^-r', s+params.dr, p2, '--r', s+params.dr, p3,'>-g', ...
    [zer zer], [0 m], '--k', ...
    ...% [zer zer]+out.LSE(ti), [0 m], ':k', ...
    LineWidth=1.5);
% plot(s, -1 + min(params.alpha2, (exp(abs(params.alpha2*(s-0.5*params.dr).^params.alpha3)))));

text(0 + params.dS/2, m, ...
    sprintf('p1\\_1: %1.2e\np2\\_1: %1.2e', out.p1_1(ti), out.p2_1(ti)));
xlim([s(1), s(end)]);
xlabel('Strain from attachment ($\mu$m)', Interpreter='latex');
ylabel('States probability', Interpreter='latex');
legend('A1', 'A2', 'A2+stroke');
if ~isequal(yl, [0, 1])
    ylim(yl);
end

if isempty(hObject)
    ylim([0 max([p1(:); p2(:); p3(:)])])
end



end

function rescalePlot(event)
    g = gca;
    ylim('auto');
end