function writeBinFile(fn, data, rcFlag, byteOrder, dataType)
% Creator John Merryman - EMI
% Date: 02 Sep 2009
% usage: writeBinFile(FN, DATA, [RCFLAG], [BYTEORDER], [DATATYPE])
%
% Write matrix to file in row-major order.
% Complex data is written interleaving I,Q
%
% FN        : (input) output file name
% DATA      : (input) data matrix (2D).
% RCFLAG    : (input) 1=real data (default); 2=complex data
% BYTEORDER : string specifying desired byte-ordering convention
%             (e.g. 'ieee-be','ieee-le'. See fopen documentation.)
% DATATYPE  : string specifying thr typr of DATA's elements 
%             (e.g. 'float32'. See fopen documentation.)
% Example:
%     writeBinFile('/home/myfile.bin', matrix, 1, 'ieee-be', 'float32');
%

if nargin < 2
    error('Filename and data matrix must be specified, exiting.');
end

if nargin < 3
  rcFlag = 1;
end

if nargin < 4
  byteOrder = 'ieee-be';
end

if nargin < 5
  dataType = 'float32';
end

[nrows,ncols] = size(data);

if rcFlag == 2
   outData = zeros(nrows,2*ncols);
   outData(:,1:2:end) = real(data);
   outData(:,2:2:end) = imag(data);
else
   outData = data;
end

fid = fopen(fn, 'w', byteOrder);
count = fwrite(fid, transpose(outData), dataType);
fclose(fid);

if count ~= size(outData,1)*size(outData,2)
    error('Problem writing data');
end
