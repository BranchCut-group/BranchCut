
%data = readBinFile('longs.152x458.phase', 152, 1, 'ieee-be', 'int8');
data = readBinFile('isola.157x458.phase', 157, 1, 'ieee-be', 'int8');
pha = pi * (data + 128) / 128 - pi; % wrapped phase in radians
