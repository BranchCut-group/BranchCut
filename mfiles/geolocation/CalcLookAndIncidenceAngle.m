function [lookAngle, incAngle] = CalcLookAndIncidenceAngle(sarParameters, ...
    sample, line, height, elliSemiMajor, elliFlat)
%
% Creator John Peter Merryman Boncori - sarmap 
% Date: 11 Jan 2017
%
% Usage: 
%
% lookAngle, incAngle = CalcLookAndincidenceAngle(sarParameters, sample, line, [height], ...
%                           [elliSemiMajor], [elliFlat]);
%
% Compute satellite look and incidence angles for a given pixel 
% (sample,line,height):
% 
%
% sarParameters     : (input)  Raw (or SLC) data parameter struct
% sample            : (input)  Raw or SLC range sample number
% line              : (input)  Raw or SLC azimuth line number
% height            : (input)  pixel height (default = 500 m)
% ellSemiMaj        : (input)  ellipsoid semi major axis (m) (default = WGS84 parameters)
% ellFlat           : (input)  ellipsoid flattening (default = (default = WGS84 parameters))
% lookAngle         : (output) Look angle (degrees)
% incAngle          : (output) Incidence angle (degrees)
%
% Example:
%     S = ParseISPpar('19960410.rslc.par');
%     [lookAngle, incAngle] = CalcLookAndIncidenceAngle(S, S.nsa/2, S.nli/2);
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

% Compute coordinates of the point on ground for which the given sample/line
% represent the zero-Doppler radar coordinates
elliEccSquared = elliFlat*(2-elliFlat);
[latP, lonP] = CalcSARLatLon(sample, line, height, sarParameters);

% Convert to Cartesian coordinates.
[xP, yP, zP] = ell2xyz(deg2rad(latP), deg2rad(lonP), height, elliSemiMajor, elliEccSquared);
PointXYZ = [xP; yP; zP];

% Compute satellite position and velocity polynomials from state vectors
azimuthTime = sarParameters.tstart + (line - 1) / sarParameters.PRF;
if ~isfield(sarParameters, 'satelliteTrajectory')
    tstop = sarParameters.tstart + (sarParameters.nli - 1) / sarParameters.PRF;
    XPoly = CalcSatTrajectory(sarParameters.stateVectors, sarParameters.tstart, tstop);
else
    XPoly = sarParameters.satelliteTrajectory.XPoly;
end
XSat = polyvalSV(XPoly, azimuthTime); 

% Line of sight vector
LOS = PointXYZ - XSat;

% Off-nadir (look) angle
cosLookAngle = dot(XSat,LOS);
sinLookAngle = norm(cross(XSat,LOS));
lookAngle = rad2deg(pi - atan2(sinLookAngle,cosLookAngle));

% Off-local vertical (incidence) angle
cosIncAngle = dot(PointXYZ,LOS);
sinIncAngle = norm(cross(PointXYZ,LOS));
incAngle = rad2deg(pi - atan2(sinIncAngle,cosIncAngle));
