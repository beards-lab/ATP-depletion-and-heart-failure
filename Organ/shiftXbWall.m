function [w, W, lost] = shiftXbWall(w, dirn, W)
%SHIFTXBWALL  Slide one wall's strain window by a third of the grid.
%
% PROVENANCE: [GLUE], reproducing what Model/evaluateModel.m does in its
% movingWindow handling: the probability arrays move by floor(ss/3) bins, the
% vacated bins are zeroed, and the pivot LXBpivot moves with them. `lost` is the
% attached mass that left the grid (zero unless heads sat in the vacated bins, #52).
%   dirn = 1: the window moves right;  2: left.

ss = W.P.ss;  Ns = W.P.NumberOfStates;  nds = W.nds;
m0 = W.P.dS*sum(w(1:W.nBins));
for k = 0:Ns-1
    b = k*ss;
    if dirn == 1
        w(b+1:b+ss-nds) = w(b+1+nds:b+ss);  w(b+ss-nds+1:b+ss) = 0;
    else
        w(b+1+nds:b+ss) = w(b+1:b+ss-nds);  w(b+1:b+nds) = 0;
    end
end
if dirn == 1, W.P.LXBpivot = W.P.LXBpivot - nds*W.P.dS*2; else, W.P.LXBpivot = W.P.LXBpivot + nds*W.P.dS*2; end
lost = m0 - W.P.dS*sum(w(1:W.nBins));
end
