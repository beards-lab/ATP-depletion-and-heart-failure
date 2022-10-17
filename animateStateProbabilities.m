function animateStateProbabilities(out, params, times)
% Animate the probabilities from the EvaluateModel output (out). Make sure
% the 'ValuesInTime' option was true to generate time output. Optionally
% input times to interpolate animation values at given timepoints.

% times2 = linspace(out.t(1), out.t(end), 1000);
if nargin < 3
    times = out.t;
end
ut = find(out.t ~= circshift(out.t, 1));
PU = interp1(out.t(ut), out.PU(ut, :), times);


figure(1);clf;
% N = opts.N;

% params.dS = opts.Slim/opts.N;
%     params.Slim = opts.Slim;
%     params.s = (0:1:opts.N)*params.dS; % strain 
%     params.s = (-para.N:opts.N)*params.dS; % strain space

for i = 1:size(PU, 1)
    
    ss = length(params.s);
    p1 = PU(i, 1:ss);
    p2 = PU(i, ss +1:2*ss);
    p3 = PU(i, 2*ss+1:3*ss);
    NR = PU(i, 3*ss + 1);
    NP = PU(i, 3*ss +2);
    SL = PU(i, 3*ss +3);    
    LSe = PU(i, 3*ss +4);    
    subplot(211);cla;hold on;
    plot(params.s, p1, params.s, p2, params.s, p3);
    plot([out.Vxb(i) out.Vxb(i)]/100, [0 max([p1, p2, p3])])
%     ylim([0 500]);
    subplot(212);cla;hold on;
    plot(out.t, out.SL)
    plot(out.t, out.LXB)
    plot([times(i) times(i)], [2 2.2])
    xlim([out.t(1), out.t(end)])
    
%     drawnow;
    pause(1/30);
end