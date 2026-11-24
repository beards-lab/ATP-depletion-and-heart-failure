function [ATP, ADP] = atpLadder()
%ATPLADDER  The shared nucleotide ladder (mM). Free ADP from the TWorld lowAtp analysis:
%   8 mM -> human CK equilibrium (getCKEquilibriumADP(8) = 0.0937), 4 -> 0.182, 2 -> 0.109
%   (lowest published ischaemic free ADP at 50 % / 25 % ATP); 6 mM interpolated (0.138).
ATP = [8 6 4 2];
ADP = [0.0937 0.138 0.182 0.109];
end
