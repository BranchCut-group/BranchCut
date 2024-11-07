
%
% File        : 20031203_20040211.cut.pha
% Description : Interferometric phase
% Size        : 2065 rows x 1959 columns
% Format      : 32 bit floating point, Big Endian byte ordering
%

%
% File        : 20031203_20040211.cut.cor
% Description : Interferometric phase coherence magnitude
% Size        : 2065 rows x 1959 columns
% Format      : 32 bit floating point, Big Endian byte ordering
%

% The contents of each file can be read into a matrix and displayed as follows:

numOfColumns = 1959;

pha = readBinFile('20031203_20040211.cut.pha', numOfColumns, 1);
figure; imagesc(pha, [-pi pi]); colormap(jet)

cor = readBinFile('20031203_20040211.cut.cor', numOfColumns, 1);
figure; imagesc(cor, [0 1]); colormap(jet)
