function data = readBinFile(fn, width, rcFlag, byteOrder, dataType, r1, rN, c1, cN)
%
% Creator John Merryman - EMI
% Date: 02 Sep 2009
%
% usage: DATA = readBinFile(FN, WIDTH, [RCFLAG], [BYTEORDER], [DATATYPE],
%                           [R1], [RN], [C1], [CN])
%
% Read Rows R1 to RN and columns C1 to CN from FN into matrix DATA.
% Row and column numbers follow matlab's indexing convention
% (row,col = 1,1 corresponds to the first element of DATA).
% Elements are placed in DATA in row-major order.
%
% FN        : (input)  input file name.
% WIDTH     : (input)  samples per line in input file.
% RCFLAG    : (input)  1=real data (default); 2=complex data (assumed IQ interleaved)
% BYTEORDER : (input)  string specifying desired byte-ordering convention
%                      (default = 'ieee-be'. See fopen documentation.).
% DATATYPE  : (input)  string specifying thr typr of DATA's elements 
%                      (default = 'float32'. See fopen documentation.).
% R1,RN     : (input)  Read rows R1 to RN included (default = all rows). 
% C1,CN     : (input)  Read columns C1 to CN included (default = all columns). 
% DATA      : (output) data matrix.
%
% Example:
%   [matrix]=readBinFile('/home/myfile.bin',100,1,'ieee-be','float32',1,50,1,50);
%
  if nargin < 2
    error('Filename and width must be specified, exiting.');
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

  if nargin < 6
      r1=1;
  end

  if nargin < 7
      rN = 0;
  end

  if nargin < 8
      c1 = 1;
  end

  if nargin < 9
      cN = width;
  end

  % Number of bytes per matrix element
  switch dataType
    case {'float32'}
      elemChars = 4;
    case {'int16'}
      elemChars = 2;
    case {'int8'}
      elemChars = 1;
    otherwise
      error ('Invalid datatype specified.');
  end

  % If elements are complex, double size
  if rcFlag == 2
     elemChars = 2*elemChars;
  end

  % If number of rows not given, estimate from file size
  if rN == 0
      finfo = dir(fn);
      rN = finfo.bytes/width/elemChars;
  end

  Ncols = cN - c1 + 1;
  Nrows = rN - r1 + 1;

  % Offset in bytes to start of read
  offset = ((r1-1) * width + c1 - 1) * elemChars;

  % Bytes to skip at each read
  skip = (width - Ncols) * elemChars;

  % Read in data
  fid = fopen(fn, 'r', byteOrder);
  fseek(fid, offset, 'bof');
  %data = zeros(Nrows,rcFlag*Ncols);
  data = fread(fid, [rcFlag*Ncols Nrows], [num2str(rcFlag * Ncols), '*', dataType], skip, byteOrder);
  
  % transpose to store in row-major format (fread stores in column-major)
  data = data.'; 

  if rcFlag == 2
     data = complex(data(:,1:2:end),data(:,2:2:end));
  end

  fclose(fid);

end
