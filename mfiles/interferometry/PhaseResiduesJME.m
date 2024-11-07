function residueMap = PhaseResiduesJME(phase)
%
% Creator John Merryman, DTU Space
% Date: 05 Mar. 2018
%
% Usage:
%
% residueMap = PhaseResiduesJME(phase)
%
% Compute a map of phase residues.
% NOTE: if input matrix size is nrows x ncols -> residue map is nrows + 1 x ncols + 1
%
% phase      : (input)  Interferometric phase matrix (real values in rad)where.
% residueMap : (output) Map with 1,-1 (positive and negative residues) and 0 (no residue).
%
% Example:
%    pha = angle(complexInterferogram);
%    residues = PhaseResiduesJME(pha); 
%


[nrows, ncols] = size(phase);

rowDiff = wrap(diff(phase, 1, 1));
colDiff = wrap(diff(phase, 1, 2));

% Clockwise circuitation
residueMap = nan(nrows + 1, ncols + 1);
residueMap(2 : end - 1, 2 : end - 1) = round((diff(rowDiff,1,2) - diff(colDiff,1,1)) / 2 / pi);
