function [C, parts, fM] = stairsCost(t, F, velocitytable, fD)
%STAIRSCOST  Block-balanced feature cost of a staircase simulation vs data.
%
%   [C, PARTS, FM] = stairsCost(T, F, VELOCITYTABLE, FD)
%
%   FD is the data feature struct from extractStairsAttributes(..., isData=true).
%   Each of the four per-step response blocks contributes
%
%       C_block = mean_k( (m_k - d_k)^2 ) / rms_k(d_k)^2
%
%   so every block is dimensionless, scaled by its own measured amplitude, and
%   equal to 1.0 when the model is wrong by one full block amplitude. This
%   deliberately gives the four PHASES of the step response equal say --
%   the raw MSE instead gives ~90 % of its weight to the holds simply because
%   they occupy 90 % of the samples, and most of that is measurement noise.
%
%   Blocks: peak (elastic/T1), nadir (phase-2 undershoot), plat (envelope),
%           rec (delayed recovery within the hold).

fM = extractStairsAttributes(t, F, velocitytable, 'st_', struct(), false);

blocks = {'st_peak','st_nadir','st_plat','st_rec'};
parts  = struct();
C      = 0;
for i = 1:numel(blocks)
    b = blocks{i};
    d = fD.(b)(:); m = fM.(b)(:);
    n = min(numel(d), numel(m));
    c = mean((m(1:n) - d(1:n)).^2) / mean(d(1:n).^2);
    parts.(b) = c;
    C = C + c;
end
end
