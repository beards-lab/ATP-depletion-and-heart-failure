% draw interactive plot in time
figure(11);clf;
makeplot([], [], out);
h = uicontrol('style','slider','units','pixel','position',[20 20 300 20]);
h.Callback = @(hObject, event) makeplot(hObject, event,out);
% addListener(h, 'ContinuousValueChange', @(hObject, event) makeplot(hObject, event,out));

function makeplot(hObject, event,out)
if isempty(hObject) 
    ti = 1;
    t = out.t(1);
else
    sval = get(hObject,'Value');
    t = sval*(out.t(end) - out.t(1)) + out.t(1);
    ti = find(out.t >= t, 1);
end

disp(hObject)
disp(event)

subplot(311);cla;hold on; 

% t = min(out.t);

Pus = 1 - out.p1_0 - out.p2_0 - out.p3_0;% PU substitute
leg = plot(out.t, Pus, out.t, out.p1_0, out.t, out.p2_0, out.t, out.p3_0, out.t, out.NR);
legend(leg); legend('Pu', 'P1', 'P2', 'P3', 'NR');
plot([t t], [0 1]);
text(t,0.5, num2str(t));

subplot(312);cla;hold on;
ss = 101; % params.ss;
s = -0.1:0.002:0.1;%params.s;
% ti = 1;
p1 = out.PU(ti, 1:ss);
p2 = out.PU(ti, ss+1:2*ss);
p3 = out.PU(ti, 2*ss+1:3*ss);

plot(s, p1, s, p2, s, p3);

end