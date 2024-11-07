function Vg = CalcVg(sarParameters, sample, line, height, elliSemiMajor, elliFlat)
%
% Creator John Peter Merryman Boncori - sarmap 
% Date: 14 Feb 2017
%
% Usage: 
%
% Vg = CalcVg(sarParameters, sample, line, [height], [elliSemiMajor], [elliFlat]);
%
% Compute ground velocity for a given pixel (sample,line,height):
% 
% sarParameters     : (input)  Raw (or SLC) data parameter struct
% sample            : (input)  Raw or SLC range sample number
% line              : (input)  Raw or SLC azimuth line number
% height            : (input)  pixel height (default = 500 m)
% ellSemiMaj        : (input)  ellipsoid semi major axis (m) (default = WGS84 parameters)
% ellFlat           : (input)  ellipsoid flattening (default = (default = WGS84 parameters))
% Vg                : (output) Ground velocity (m / s)
%
% Example:
%     S = ParseISPpar('19960410.rslc.par');
%     Vg = CalcVg(S, S.nsa/2, S.nli/2);
%

if nargin < 3
  error('Insufficient number of input parameters, exiting.');
end

if nargin < 4
  height  = 500; 
end

if nargin < 5
  WGS84_MAJ  = 6378137; 
  elliSemiMajor = WGS84_MAJ; 
end

if nargin < 6
  WGS84_FLAT = 1/298.257223563;
  elliFlat = WGS84_FLAT; 
end   

% Compute coordinates for two azimuth lines at the given slant range
elliEccSquared = elliFlat*(2-elliFlat);

[latP1, lonP1] = CalcSARLatLon(sample, line + 1, height, sarParameters);
[latP2, lonP2] = CalcSARLatLon(sample, line - 1, height, sarParameters);

% Convert to Cartesian coordinates.
[xP1, yP1, zP1] = ell2xyz(deg2rad(latP1), deg2rad(lonP1), height, elliSemiMajor, elliEccSquared);
Point1XYZ = [xP1; yP1; zP1];

[xP2, yP2, zP2] = ell2xyz(deg2rad(latP2), deg2rad(lonP2), height, elliSemiMajor, elliEccSquared);
Point2XYZ = [xP2; yP2; zP2];

% Distance (great circle)
Rlocal = elliSemiMajor * (1 - elliFlat * sin(deg2rad(latP1))^2);
d = Rlocal * atan2(norm(cross(Point1XYZ,Point2XYZ)), dot(Point1XYZ,Point2XYZ));

% Ground (beam) velocity
Vg = d / 2 * sarParameters.PRF;
